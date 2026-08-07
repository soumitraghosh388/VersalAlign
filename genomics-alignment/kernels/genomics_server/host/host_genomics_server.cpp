// genomics_server host — VCK5000  (RESIDENT server model)
//
// Difference vs genomics_dma (per-launch): the two server CUs are launched
// ONCE and each processes NUM_ITERS batches via a DDR doorbell. Per batch
// the host pays only a doorbell write + completion poll — NOT a kernel-launch
// (ERT/XGQ scheduler) round-trip. This isolates and removes the per-batch
// launch overhead measured in the per-launch model.
//
//   dist_server(batch_in, ctrl_in, num_iters, 8 row streams)
//       loop b: wait until ctrl_in[0] >= b+1, then stream batch b in.
//   collect_server(scores_out, ctrl_out, num_iters, 50 col streams)
//       loop b: drain one batch, store it, set ctrl_out[0] = b+1.
//
// Data layout (from gen_multi_data.py, same as genomics_dma_swg):
//   input_r.txt   — NUM_ITERS x BUF_SIZE int32 words for row r
//   golden_c.txt  — per column c, NUM_ITERS x [count, scores...] in
//                   bottom-row-first order, refs strided by NUM_COLS.
//
// RING=2 double-buffer pipeline:
//   slot for batch b  = (b-1) % RING
//   While HW processes batch b from slot (b-1)%RING, the host fills
//   and H2Ds batch b+1 into slot b%RING — PCIe DMA and PL/AIE overlap.
//   dist_server picks up doorbell b+1 the moment it finishes batch b.
//
// Usage: ./genomics_server.exe <xclbin> <data_dir> [bdf]

#include <cstring>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <chrono>

#include "xrt.h"
#include "experimental/xrt_kernel.h"

#define NUM_ROWS    8
#define NUM_COLS    50
#define BUF_SIZE    2048    // int32 per row per batch
#define SCORE_BUF   64     // int32 per column output slot
#define NUM_ITERS   1000
#define TOTAL_ITER  1   // for calculating power, increase the value to 100, otherwise 1
#define READ_LEN    128
#define REF_LEN     128
#define REF_PER_ROW 200
#define RING        2       // must match pl/gen_servers.py RING

using clk = std::chrono::steady_clock;
static double ms_since(clk::time_point t0) {
    return std::chrono::duration<double, std::milli>(clk::now() - t0).count();
}

static std::vector<int32_t> read_data_file(const std::string& path) {
    std::vector<int32_t> data;
    std::ifstream f(path);
    if (!f.is_open()) { std::cerr << "Cannot open " << path << std::endl; return data; }
    std::string line;
    while (std::getline(f, line))
        if (!line.empty()) data.push_back(std::stoi(line));
    return data;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <XCLBIN> <data_dir> [bdf]" << std::endl;
        return 1;
    }
    std::string xclbin_file = argv[1];
    std::string data_dir    = argv[2];
    std::string bdf = (argc > 3) ? std::string(argv[3]) : std::string("0000:82:00.1");

    std::cout << "genomics_server " << NUM_ROWS << "x" << NUM_COLS
              << " = " << NUM_ROWS * NUM_COLS << " tiles, RESIDENT"
              << "  num_iters=" << NUM_ITERS << std::endl;

    // ----- load input files (input_r.txt, one per row) -----
    std::vector<std::vector<int32_t>> inputs(NUM_ROWS);
    for (int r = 0; r < NUM_ROWS; r++) {
        inputs[r] = read_data_file(data_dir + "/input_" + std::to_string(r) + ".txt");
        if (inputs[r].empty()) {
            std::cerr << "Failed to load input_" << r << ".txt" << std::endl;
            return 1;
        }
        while ((int)inputs[r].size() < NUM_ITERS * BUF_SIZE)
            inputs[r].push_back(0);
    }

    // ----- load golden files (golden_c.txt, one per column) -----
    std::vector<std::vector<int32_t>> goldens(NUM_COLS);
    for (int c = 0; c < NUM_COLS; c++) {
        goldens[c] = read_data_file(data_dir + "/golden_" + std::to_string(c) + ".txt");
        if (goldens[c].empty()) {
            std::cerr << "Failed to load golden_" << c << ".txt" << std::endl;
            return 1;
        }
    }

    // Precompute per-column, per-iteration score start offset in goldens[c].
    // goldens[c] layout: count_s  score_0 ... score_{count-1}  (repeated NUM_ITERS times).
    // gold_off[c][s] = index of first score for iteration s in goldens[c].
    std::vector<std::vector<int>> gold_off(NUM_COLS, std::vector<int>(NUM_ITERS, 0));
    for (int c = 0; c < NUM_COLS; c++) {
        int pos = 0;
        for (int s = 0; s < NUM_ITERS; s++) {
            gold_off[c][s] = pos + 1;
            pos += 1 + goldens[c][pos];
        }
    }

    // ----- device + kernels -----
    auto device = xrt::device(bdf);
    std::cout << "Device: " << device.get_info<xrt::info::device::name>()
              << " (" << bdf << ")" << std::endl;
    auto uuid = device.load_xclbin(xclbin_file);

    auto dist    = xrt::kernel(device, uuid, "dist_server:{dist_server_0}");
    auto collect = xrt::kernel(device, uuid, "collect_server:{collect_server_0}");

    // ----- buffers (RING slots each) -----
    static constexpr size_t BATCH_SLOT_BYTES  = NUM_ROWS * BUF_SIZE  * sizeof(int32_t);
    static constexpr size_t SCORES_SLOT_BYTES = NUM_COLS * SCORE_BUF * sizeof(int32_t);

    auto batch_bo  = xrt::bo(device, RING * BATCH_SLOT_BYTES,  dist.group_id(0));
    auto ctrl_in   = xrt::bo(device, 16   * sizeof(int32_t),   dist.group_id(1));
    auto scores_bo = xrt::bo(device, RING * SCORES_SLOT_BYTES, collect.group_id(0));
    auto ctrl_out  = xrt::bo(device, 16   * sizeof(int32_t),   collect.group_id(1));

    auto* ctrl_in_map  = ctrl_in.map<int32_t*>();
    auto* ctrl_out_map = ctrl_out.map<int32_t*>();
    ctrl_in_map[0]  = 0;   ctrl_in.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    ctrl_out_map[0] = 0;   ctrl_out.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    std::vector<int32_t> zeros(RING * NUM_COLS * SCORE_BUF, 0);
    scores_bo.write(zeros.data());
    scores_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    std::vector<int32_t> batch_slot(NUM_ROWS * BUF_SIZE, 0);
    std::vector<int32_t> scores(NUM_COLS * SCORE_BUF, 0);

    double total_ms;
    std::vector<std::vector<int32_t>> all_scores;
    std::vector<double> lat;

    for(int pow_count=0; pow_count<TOTAL_ITER;pow_count++){

        // Reset doorbells so re-launched servers start from a clean state.
        ctrl_in_map[0]  = 0;   ctrl_in.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        ctrl_out_map[0] = 0;   ctrl_out.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        scores_bo.write(zeros.data());
        scores_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE);

        // ----- launch the two servers ONCE (resident for NUM_ITERS batches) -----
        std::cout << "Launching resident servers (1 launch each)..." << std::endl;
        auto run_collect = collect(scores_bo, ctrl_out, NUM_ITERS);
        auto run_dist    = dist(batch_bo, ctrl_in, NUM_ITERS);

        // ----- RING=2 double-buffer doorbell loop -----
        // b is 1-indexed (doorbell protocol); iteration data index is b-1 (0-based).
        lat.clear();
        lat.reserve(NUM_ITERS);
        // Store scores for all iterations; validated after timing window closes.
        all_scores.assign(NUM_ITERS, std::vector<int32_t>(NUM_COLS * SCORE_BUF));
        auto t_all = clk::now();

        // Prime: fill and H2D iteration 0 into slot 0, ring doorbell 1.
        for (int r = 0; r < NUM_ROWS; r++)
            std::copy(inputs[r].begin(), inputs[r].begin() + BUF_SIZE,
                    batch_slot.begin() + r * BUF_SIZE);
        batch_bo.write(batch_slot.data(), BATCH_SLOT_BYTES, 0);
        batch_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE, BATCH_SLOT_BYTES, 0);
        ctrl_in_map[0] = 1;
        ctrl_in.sync(XCL_BO_SYNC_BO_TO_DEVICE);

        for (int b = 1; b <= NUM_ITERS; b++) {
            auto t0 = clk::now();

            // Step A: H2D iteration b into next slot while HW processes iteration b-1.
            if (b < NUM_ITERS) {
                int next_slot = b % RING;
                int src_off   = b * BUF_SIZE;
                for (int r = 0; r < NUM_ROWS; r++)
                    std::copy(inputs[r].begin() + src_off,
                            inputs[r].begin() + src_off + BUF_SIZE,
                            batch_slot.begin() + r * BUF_SIZE);
                batch_bo.write(batch_slot.data(), BATCH_SLOT_BYTES,
                            (size_t)next_slot * BATCH_SLOT_BYTES);
                batch_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE, BATCH_SLOT_BYTES,
                            (size_t)next_slot * BATCH_SLOT_BYTES);
                ctrl_in_map[0] = b + 1;
                ctrl_in.sync(XCL_BO_SYNC_BO_TO_DEVICE);
            }

            // Step B: Poll for batch b completion.
            do {
                ctrl_out.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
            } while (ctrl_out_map[0] < b);

            // Step C: D2H only the slot that batch b wrote.
            int curr_slot = (b - 1) % RING;
            scores_bo.sync(XCL_BO_SYNC_BO_FROM_DEVICE, SCORES_SLOT_BYTES,
                        (size_t)curr_slot * SCORES_SLOT_BYTES);
            scores_bo.read(all_scores[b - 1].data(), SCORES_SLOT_BYTES,
                        (size_t)curr_slot * SCORES_SLOT_BYTES);

            lat.push_back(ms_since(t0));
        }

        total_ms = ms_since(t_all);
        run_dist.wait();
        run_collect.wait();
    }

    // ----- validation (outside timing window) -----
    int total_fail = 0;
    for (int s = 0; s < NUM_ITERS; s++) {
        const std::vector<int32_t>& scores = all_scores[s];
        int iter_fail = 0;
        for (int c = 0; c < NUM_COLS; c++) {
            const int32_t* col     = &scores[c * SCORE_BUF];
            int score_count        = col[0];
            int expected_count     = goldens[c][gold_off[c][s] - 1];
            if (score_count != expected_count) {
                iter_fail++;
                if (total_fail + iter_fail <= 10)
                    std::cerr << "  iter=" << s << " col=" << c
                              << ": count=" << score_count
                              << " expected=" << expected_count << "\n";
            } else {
                for (int i = 0; i < expected_count; i++) {
                    int got = col[1 + i];
                    int exp = goldens[c][gold_off[c][s] + i];
                    if (got != exp) {
                        iter_fail++;
                        if (total_fail + iter_fail <= 10)
                            std::cerr << "  iter=" << s << " col=" << c
                                      << " score[" << i << "]: got=" << got
                                      << " exp=" << exp << "\n";
                    }
                }
            }
        }
        total_fail += iter_fail;
        if (s == 0 || iter_fail)
            std::cout << "  iter " << s << ": "
                      << (iter_fail ? "FAIL" : "ok") << "  "
                      << std::fixed << std::setprecision(3) << lat[s] << " ms"
                      << std::endl;
    }

    // ----- metrics -----
    std::sort(lat.begin(), lat.end());
    double sum = 0; for (double x : lat) sum += x;
    double mn = lat.front(), md = lat[lat.size() / 2], mx = lat.back();
    double ss = (lat.size() > 1) ? (sum - lat[0]) / (lat.size() - 1) : lat[0];
    double gcups = (double)NUM_ITERS * NUM_ROWS * REF_PER_ROW * READ_LEN * REF_LEN
                   / (total_ms * 1e6);

    std::cout << "\n=== Results (resident, num_iters=" << NUM_ITERS << ") ===" << std::endl;
    std::cout << "Tiles:            " << NUM_ROWS * NUM_COLS << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Per-iter (doorbell): min " << mn << " / median " << md
              << " / max " << mx << " ms" << std::endl;
    std::cout << "Steady-state mean (excl. iter 0): " << ss << " ms" << std::endl;
    std::cout << "Total wall time:  " << total_ms << " ms for "
              << NUM_ITERS << " iters" << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Throughput:       " << gcups << " GCUPS  [wall clock]" << std::endl;
    std::cout << (total_fail == 0 ? "\nPASS\n" : "\nFAIL\n") << std::endl;

    return total_fail == 0 ? 0 : 1;
}
