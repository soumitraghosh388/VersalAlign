#include <iostream>
#include <fstream>
#include <string>
#include "graph.hpp"

bitapGraph ba_graph;

#if defined(__AIESIM__) || defined(__X86SIM__)
int main(int argc, char** argv) {
    ba_graph.init();
    ba_graph.run(1);
    ba_graph.end();

    std::ifstream golden_file, aie_file;
    golden_file.open("../data/out_score.txt");
    if (golden_file.fail()) {
        std::cerr << "Error opening golden file" << std::endl;
        return -1;
    }

#if defined(__X86SIM__)
    aie_file.open("x86simulator_output/output_score.txt");
#else
    aie_file.open("aiesimulator_output/output_score.txt");
#endif
    if (aie_file.fail()) {
        std::cerr << "Error opening AIE output file" << std::endl;
        return -1;
    }

    std::string line_golden, line_aie;
    int match_count = 0, mismatch_count = 0, idx = 0;
    while (getline(golden_file, line_golden)) {
        if (!getline(aie_file, line_aie)) {
            std::cerr << "AIE output too short at idx " << idx << std::endl;
            mismatch_count++;
            break;
        }
        while (!line_aie.empty() && line_aie[0] == 'T')
            getline(aie_file, line_aie);

        int golden_val = std::stoi(line_golden);
        int aie_val = std::stoi(line_aie);
        if (golden_val != aie_val) {
            std::cerr << "MISMATCH idx " << idx << ": expected "
                      << golden_val << " got " << aie_val << std::endl;
            mismatch_count++;
        } else {
            match_count++;
        }
        idx++;
    }

    if (mismatch_count == 0)
        std::cout << "PASS: Bitap - " << match_count << " scores matched" << std::endl;
    else
        std::cout << "FAIL: Bitap - " << mismatch_count << " mismatches" << std::endl;

    return mismatch_count;
}
#endif
