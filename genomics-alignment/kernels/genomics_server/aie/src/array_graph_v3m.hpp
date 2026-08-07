#ifndef ARRAY_GRAPH_V3M_HPP
#define ARRAY_GRAPH_V3M_HPP

#include <adf.h>
using namespace adf;

#include "tile_bottom_v3.hh"
#include "tile_upper_v3.hh"
#include "kernel_config.h"

class ArrayGraphV3M : public graph {
private:
    kernel tiles[NUM_ROWS][NUM_COLS];

public:
    input_plio  plio_in[NUM_ROWS];
    output_plio plio_out[NUM_COLS];

    ArrayGraphV3M() {
        for (int r = NUM_ROWS - 1; r >= 0; r--) {
            std::string rid = std::to_string(r);

            plio_in[r] = input_plio::create("in_" + rid, plio_32_bits,
                             "data/input_" + rid + ".txt");

            for (int c = 0; c < NUM_COLS; c++) {
                std::string cid = std::to_string(c);

                if (r == NUM_ROWS - 1) {
                    tiles[r][c] = kernel::create_object<TileBottomV3>(r, c);
                    source(tiles[r][c]) = "tile_bottom_v3.cc";
                    dimensions(tiles[r][c].in[0]) = {BUF_SIZE};
                    dimensions(tiles[r][c].out[0]) = {SCORE_BUF};
                } else {
                    tiles[r][c] = kernel::create_object<TileUpperV3>(r, c);
                    source(tiles[r][c]) = "tile_upper_v3.cc";
                    dimensions(tiles[r][c].in[0]) = {BUF_SIZE};
                    dimensions(tiles[r][c].in[1]) = {SCORE_BUF};
                    dimensions(tiles[r][c].out[0]) = {SCORE_BUF};
                    connect(tiles[r + 1][c].out[0], tiles[r][c].in[1]);
                }

                runtime<ratio>(tiles[r][c]) = 0.9;
                location<kernel>(tiles[r][c]) = tile(c, r);
                connect(plio_in[r].out[0], tiles[r][c].in[0]);

                if (r == 0) {
                    plio_out[c] = output_plio::create("out_" + cid, plio_32_bits,
                                      "output_col_" + cid + ".txt");
                    connect(tiles[0][c].out[0], plio_out[c].in[0]);
                }
            }
        }
    }
};

#endif
