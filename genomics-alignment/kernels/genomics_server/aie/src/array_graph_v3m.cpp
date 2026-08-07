#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "array_graph_v3m.hpp"

ArrayGraphV3M array_graph_v3m;

#if defined(__AIESIM__) || defined(__X86SIM__)

// Read next integer from file, skipping aiesimulator 'T ...' timestamp lines.
static int read_next_int(std::ifstream& f) {
    std::string line;
    while (getline(f, line)) {
        if (!line.empty() && line[0] != 'T')
            return std::stoi(line);
    }
    return -1;
}

int main(int argc, char** argv) {
    array_graph_v3m.init();
    array_graph_v3m.run(NUM_ITERS);
    array_graph_v3m.end();

    // Load per-column goldens: golden_c.txt for c = 0 .. NUM_COLS-1.
    // Format: count \n scores in bottom-row-first order, refs strided by NUM_COLS.
    // This matches write_column_goldens() in gen_multi_data.py.
    std::vector<std::vector<int>> golden(NUM_COLS);
    for (int c = 0; c < NUM_COLS; c++) {
        std::ifstream gf("../data/golden_" + std::to_string(c) + ".txt");
        if (gf.fail()) {
            std::cerr << "Error opening golden_" << c << ".txt" << std::endl;
            return -1;
        }
        int count; gf >> count;
        golden[c].resize(count);
        for (int i = 0; i < count; i++) gf >> golden[c][i];
    }

    int match_count = 0, mismatch_count = 0;

    // Verify all NUM_COLS column output files.
    // Column c output order (TileUpperV3 bottom-to-top aggregation):
    //   [total_count,
    //    row(ROWS-1) scores for refs c, c+COLS, c+2*COLS, ...,
    //    row(ROWS-2) scores for refs c, c+COLS, ...,
    //    ...
    //    row0 scores for refs c, c+COLS, ...]
    // This is the same order written by gen_multi_data.py write_column_goldens().
    for (int c = 0; c < NUM_COLS; c++) {
        std::string sim_path;
#if defined(__X86SIM__)
        sim_path = "x86simulator_output/output_col_" + std::to_string(c) + ".txt";
#else
        sim_path = "aiesimulator_output/output_col_" + std::to_string(c) + ".txt";
#endif
        std::ifstream af(sim_path);
        if (af.fail()) {
            std::cerr << "Error opening " << sim_path << std::endl;
            return -1;
        }

        int sim_count  = read_next_int(af);
        int gold_count = (int)golden[c].size();

        if (sim_count != gold_count) {
            std::cerr << "COUNT MISMATCH col=" << c << ": sim=" << sim_count
                      << " expected=" << gold_count << std::endl;
            mismatch_count += gold_count;
            for (int i = 0; i < sim_count; i++) read_next_int(af);
            continue;
        }

        int refs_per_col = gold_count / NUM_ROWS;  // refs handled by one row for this col
        for (int i = 0; i < gold_count; i++) {
            int aie_val = read_next_int(af);
            int exp_val = golden[c][i];
            if (aie_val != exp_val) {
                if (mismatch_count < 10) {
                    int row_from_bottom = i / refs_per_col;
                    int row = NUM_ROWS - 1 - row_from_bottom;
                    int ref = c + (i % refs_per_col) * NUM_COLS;
                    std::cerr << "MISMATCH col=" << c << " row=" << row
                              << " ref=" << ref << ": expected " << exp_val
                              << " got " << aie_val << std::endl;
                }
                mismatch_count++;
            } else {
                match_count++;
            }
        }
    }

    if (mismatch_count == 0)
        std::cout << "PASS: " << match_count << " scores matched" << std::endl;
    else
        std::cout << "FAIL: " << mismatch_count << " mismatches" << std::endl;
    std::cout << "ArrayV3 " << NUM_ROWS << "x" << NUM_COLS << " completed." << std::endl;

    return mismatch_count;
}
#endif
