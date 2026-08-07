// SWG Scalar Host Program — VCK5000 Hardware Test
//
// Reads pre-packed PLIO data files (seq1.txt, seq2.txt, out_score_sw.txt),
// sends sequences to AIE via PL mm2s kernels, reads alignment scores from
// PL s2mm, and validates against golden reference.
//
// Usage:
//   ./swg_scalar <xclbin> <data_dir> [debug]
//
// Example:
//   ./swg_scalar build.hw/swg_scalar.xclbin aie/data
//   ./swg_scalar build.hw/swg_scalar.xclbin ../../test_data/real_phix

#include <cstring>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <chrono>

#include "xrt.h"
#include "experimental/xrt_kernel.h"

// Read integers from a text file (one per line)
std::vector<int32_t> read_data_file(const std::string& path) {
    std::vector<int32_t> data;
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "ERROR: Cannot open " << path << std::endl;
        return data;
    }
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty()) {
            data.push_back(std::stoi(line));
        }
    }
    return data;
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4) {
        std::cout << "Usage: " << argv[0] << " <XCLBIN> <data_dir> [debug]" << std::endl;
        std::cout << "  data_dir should contain seq1.txt, seq2.txt, out_score_sw.txt" << std::endl;
        return EXIT_FAILURE;
    }

    std::string xclbin_file = argv[1];
    std::string data_dir = argv[2];
    bool debug = (argc == 4 && std::string(argv[3]) == "debug");

    // ---- Load data files ----
    std::cout << "Loading data from " << data_dir << "/" << std::endl;

    auto seq1_data = read_data_file(data_dir + "/seq1.txt");
    auto seq2_data = read_data_file(data_dir + "/seq2.txt");
    auto golden = read_data_file(data_dir + "/out_score_sw.txt");

    if (seq1_data.empty() || seq2_data.empty() || golden.empty()) {
        std::cerr << "ERROR: Failed to load data files" << std::endl;
        return EXIT_FAILURE;
    }

    int pair_count = seq1_data[0];
    std::cout << "  seq1: " << seq1_data.size() << " words (pair_count=" << pair_count << ")" << std::endl;
    std::cout << "  seq2: " << seq2_data.size() << " words" << std::endl;
    std::cout << "  golden: " << golden.size() << " scores" << std::endl;

    if ((int)golden.size() != pair_count) {
        std::cerr << "WARNING: golden has " << golden.size() << " scores but pair_count=" << pair_count << std::endl;
    }

    // The output buffer needs space for pair_count scores
    // Add extra space just in case
    int output_size = pair_count + 16;

    // ---- Setup device ----
    std::cout << "Open device" << std::endl;
    auto device = xrt::device(0);
    std::cout << "  name: " << device.get_info<xrt::info::device::name>() << std::endl;

    std::cout << "Load xclbin: " << xclbin_file << std::endl;
    auto uuid = device.load_xclbin(xclbin_file);

    // ---- Get PL kernel references ----
    std::cout << "Get PL kernel references" << std::endl;
    auto pl_mm2s_1 = xrt::kernel(device, uuid, "pl_mm2s:{pl_mm2s_1}");
    auto pl_mm2s_2 = xrt::kernel(device, uuid, "pl_mm2s:{pl_mm2s_2}");
    auto pl_s2mm_1 = xrt::kernel(device, uuid, "pl_s2mm:{pl_s2mm_1}");

    // ---- Allocate device buffers ----
    std::cout << "Allocate buffers" << std::endl;
    size_t seq1_bytes = seq1_data.size() * sizeof(int32_t);
    size_t seq2_bytes = seq2_data.size() * sizeof(int32_t);
    size_t out_bytes  = output_size * sizeof(int32_t);

    auto seq1_bo = xrt::bo(device, seq1_bytes, pl_mm2s_1.group_id(0));
    auto seq2_bo = xrt::bo(device, seq2_bytes, pl_mm2s_2.group_id(0));
    auto out_bo  = xrt::bo(device, out_bytes,  pl_s2mm_1.group_id(0));

    // ---- Write input data to buffers ----
    std::cout << "Write input data to buffers" << std::endl;
    seq1_bo.write(seq1_data.data());
    seq2_bo.write(seq2_data.data());

    // Zero the output buffer
    std::vector<int32_t> out_zeros(output_size, 0);
    out_bo.write(out_zeros.data());

    // ---- Sync to device ----
    std::cout << "Sync buffers to device" << std::endl;
    seq1_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    seq2_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    out_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    // ---- Run kernels ----
    std::cout << "Run kernels (" << pair_count << " alignment pairs)" << std::endl;
    auto t_start = std::chrono::high_resolution_clock::now();

    // Start s2mm first (consumer must be ready before producer)
    auto run_s2mm = pl_s2mm_1(out_bo, nullptr, pair_count);
    // Start both mm2s producers
    auto run_mm2s_1 = pl_mm2s_1(seq1_bo, nullptr, (int)seq1_data.size());
    auto run_mm2s_2 = pl_mm2s_2(seq2_bo, nullptr, (int)seq2_data.size());

    // Wait for completion
    run_mm2s_1.wait();
    std::cout << "  pl_mm2s_1 (seq1) completed" << std::endl;
    run_mm2s_2.wait();
    std::cout << "  pl_mm2s_2 (seq2) completed" << std::endl;
    run_s2mm.wait();
    std::cout << "  pl_s2mm_1 (scores) completed" << std::endl;

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    // ---- Read results ----
    std::cout << "Sync output buffer from device" << std::endl;
    out_bo.sync(XCL_BO_SYNC_BO_FROM_DEVICE);

    std::vector<int32_t> scores(output_size);
    out_bo.read(scores.data());

    // ---- Validate ----
    std::cout << std::endl << "=== Validation ===" << std::endl;
    int match_count = 0;
    int error_count = 0;
    int num_to_check = std::min(pair_count, (int)golden.size());

    for (int i = 0; i < num_to_check; i++) {
        if (scores[i] == golden[i]) {
            match_count++;
        } else {
            error_count++;
            if (error_count <= 10) {
                std::cout << "  MISMATCH pair " << i << ": got " << scores[i]
                          << " expected " << golden[i] << std::endl;
            }
        }
        if (debug) {
            std::cout << "  Pair " << std::setw(4) << i
                      << ": hw=" << std::setw(6) << scores[i]
                      << " golden=" << std::setw(6) << golden[i]
                      << " " << (scores[i] == golden[i] ? "OK" : "FAIL") << std::endl;
        }
    }

    // ---- Report ----
    std::cout << std::endl << "=== Results ===" << std::endl;
    std::cout << "Pairs:     " << num_to_check << std::endl;
    std::cout << "Matched:   " << match_count << std::endl;
    std::cout << "Errors:    " << error_count << std::endl;
    std::cout << "Time:      " << std::fixed << std::setprecision(2) << elapsed_ms << " ms" << std::endl;

    if (num_to_check > 0 && elapsed_ms > 0) {
        // Estimate average cells per pair (assume 100x120 for PhiX data)
        // In a real host, we'd compute this from the lens words
        double pairs_per_sec = num_to_check / (elapsed_ms / 1000.0);
        std::cout << "Throughput: " << std::fixed << std::setprecision(1)
                  << pairs_per_sec << " pairs/sec" << std::endl;
    }

    if (error_count == 0) {
        std::cout << std::endl << "PASS: SWG scalar - " << match_count
                  << " pairs matched on VCK5000 hardware" << std::endl;
    } else {
        std::cout << std::endl << "FAIL: SWG scalar - " << error_count
                  << " mismatches out of " << num_to_check << " pairs" << std::endl;
    }

    return error_count;
}
