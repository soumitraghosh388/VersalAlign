// SWG Inter-Sequence SIMD Host — VCK5000 Hardware Test
//
// Sends batched 16-pair data to AIE kernel, reads packed int16 scores.
//
// Usage:
//   ./swg_interseq <xclbin> <data_dir> [debug]

#include <cstring>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <chrono>

#include "xrt.h"
#include "experimental/xrt_kernel.h"

std::vector<int32_t> read_data_file(const std::string& path) {
    std::vector<int32_t> data;
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "ERROR: Cannot open " << path << std::endl;
        return data;
    }
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty()) data.push_back(std::stoi(line));
    }
    return data;
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4) {
        std::cout << "Usage: " << argv[0] << " <XCLBIN> <data_dir> [debug]" << std::endl;
        return EXIT_FAILURE;
    }

    std::string xclbin_file = argv[1];
    std::string data_dir = argv[2];
    bool debug = (argc == 4 && std::string(argv[3]) == "debug");

    // Load data
    auto seq1_data = read_data_file(data_dir + "/seq1.txt");
    auto seq2_data = read_data_file(data_dir + "/seq2.txt");
    auto golden_packed = read_data_file(data_dir + "/out_score.txt");

    if (seq1_data.empty() || seq2_data.empty() || golden_packed.empty()) {
        std::cerr << "ERROR: Failed to load data files" << std::endl;
        return EXIT_FAILURE;
    }

    int batch_count = seq1_data[0];
    int total_pairs = batch_count * 16;
    // Output: 8 packed int32 per batch (16 int16 scores packed 2 per int32)
    int output_words = batch_count * 8 + 16;  // extra padding

    std::cout << "Batches: " << batch_count << " (" << total_pairs << " pairs)" << std::endl;
    std::cout << "seq1: " << seq1_data.size() << " words, seq2: " << seq2_data.size() << " words" << std::endl;
    std::cout << "Golden: " << golden_packed.size() << " packed words" << std::endl;

    // Setup device
    std::cout << "Open device" << std::endl;
    auto device = xrt::device(0);
    std::cout << "Load xclbin: " << xclbin_file << std::endl;
    auto uuid = device.load_xclbin(xclbin_file);

    auto pl_mm2s_1 = xrt::kernel(device, uuid, "pl_mm2s:{pl_mm2s_1}");
    auto pl_mm2s_2 = xrt::kernel(device, uuid, "pl_mm2s:{pl_mm2s_2}");
    auto pl_s2mm_1 = xrt::kernel(device, uuid, "pl_s2mm:{pl_s2mm_1}");

    // Allocate buffers
    size_t seq1_bytes = seq1_data.size() * sizeof(int32_t);
    size_t seq2_bytes = seq2_data.size() * sizeof(int32_t);
    size_t out_bytes = output_words * sizeof(int32_t);

    auto seq1_bo = xrt::bo(device, seq1_bytes, pl_mm2s_1.group_id(0));
    auto seq2_bo = xrt::bo(device, seq2_bytes, pl_mm2s_2.group_id(0));
    auto out_bo  = xrt::bo(device, out_bytes,  pl_s2mm_1.group_id(0));

    seq1_bo.write(seq1_data.data());
    seq2_bo.write(seq2_data.data());
    std::vector<int32_t> out_zeros(output_words, 0);
    out_bo.write(out_zeros.data());

    seq1_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    seq2_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    out_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    // Run
    std::cout << "Run kernels" << std::endl;
    auto t_start = std::chrono::high_resolution_clock::now();

    auto run_s2mm = pl_s2mm_1(out_bo, nullptr, batch_count * 8);  // 8 int32 per batch
    auto run_mm2s_1 = pl_mm2s_1(seq1_bo, nullptr, (int)seq1_data.size());
    auto run_mm2s_2 = pl_mm2s_2(seq2_bo, nullptr, (int)seq2_data.size());

    run_mm2s_1.wait();
    run_mm2s_2.wait();
    run_s2mm.wait();

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    // Read results
    out_bo.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    std::vector<int32_t> results(output_words);
    out_bo.read(results.data());

    // Validate — compare packed int32 scores against golden
    int match_count = 0;
    int error_count = 0;
    int num_check = std::min((int)golden_packed.size(), batch_count * 8);

    // Unpack all scores for display
    std::vector<int16_t> hw_scores;
    std::vector<int16_t> golden_scores;

    for (int i = 0; i < num_check; i++) {
        int32_t hw_packed = results[i];
        int32_t gd_packed = golden_packed[i];

        int16_t hw_s0 = (int16_t)(hw_packed & 0xFFFF);
        int16_t hw_s1 = (int16_t)((hw_packed >> 16) & 0xFFFF);
        int16_t gd_s0 = (int16_t)(gd_packed & 0xFFFF);
        int16_t gd_s1 = (int16_t)((gd_packed >> 16) & 0xFFFF);

        hw_scores.push_back(hw_s0);
        hw_scores.push_back(hw_s1);
        golden_scores.push_back(gd_s0);
        golden_scores.push_back(gd_s1);

        if (hw_packed == gd_packed) {
            match_count += 2;
        } else {
            if (hw_s0 != gd_s0) {
                error_count++;
                if (error_count <= 10)
                    std::cout << "  MISMATCH pair " << (i * 2) << ": hw=" << hw_s0 << " golden=" << gd_s0 << std::endl;
            } else {
                match_count++;
            }
            if (hw_s1 != gd_s1) {
                error_count++;
                if (error_count <= 10)
                    std::cout << "  MISMATCH pair " << (i * 2 + 1) << ": hw=" << hw_s1 << " golden=" << gd_s1 << std::endl;
            } else {
                match_count++;
            }
        }
    }

    if (debug) {
        for (size_t i = 0; i < hw_scores.size() && i < golden_scores.size(); i++) {
            std::cout << "  Pair " << std::setw(3) << i
                      << ": hw=" << std::setw(6) << hw_scores[i]
                      << " golden=" << std::setw(6) << golden_scores[i]
                      << " " << (hw_scores[i] == golden_scores[i] ? "OK" : "FAIL") << std::endl;
        }
    }

    std::cout << std::endl << "=== Results ===" << std::endl;
    std::cout << "Pairs:     " << total_pairs << std::endl;
    std::cout << "Matched:   " << match_count << std::endl;
    std::cout << "Errors:    " << error_count << std::endl;
    std::cout << "Time:      " << std::fixed << std::setprecision(2) << elapsed_ms << " ms" << std::endl;

    if (total_pairs > 0 && elapsed_ms > 0) {
        double pairs_per_sec = total_pairs / (elapsed_ms / 1000.0);
        std::cout << "Throughput: " << std::fixed << std::setprecision(1) << pairs_per_sec << " pairs/sec" << std::endl;
    }

    if (error_count == 0)
        std::cout << std::endl << "PASS: SWG interseq - " << match_count << " pairs matched on VCK5000" << std::endl;
    else
        std::cout << std::endl << "FAIL: SWG interseq - " << error_count << " mismatches" << std::endl;

    return error_count;
}
