// v3 Bottom tile: compute independently, then output scores.
// Same as v1 — bottom tile doesn't wait for anyone.

#include "tile_bottom_v3.hh"
#include <aie_api/aie.hpp>
#include "kernel_config.h"
#include "compute/compute.hh"

TileBottomV3::TileBottomV3(int row, int col) : my_row(row), my_col(col) {}

void TileBottomV3::run(
    adf::input_async_buffer<int32_t>& __restrict batch_in,
    adf::output_async_buffer<int32_t>& __restrict scores_out
) {
    // Phase 1: Compute (only needs batch_in)
    batch_in.acquire();
    auto in_ptr = (int32_t*)batch_in.data();
    int32_t total_items = in_ptr[0];

    SwProfile prof;
    compute_init(in_ptr, total_items, prof);

    // Compute all my scores into local array
    int32_t my_scores[64];  // max refs per tile
    int my_count = 0;
    for (int i = my_col; i < total_items; i += NUM_COLS) {
        my_scores[my_count] = compute_score(in_ptr, total_items, i, my_row, my_col, prof);
        my_count++;
    }
    batch_in.release();

    // Phase 2: Output scores (after compute is fully done)
    scores_out.acquire();
    auto out_ptr = (int32_t*)scores_out.data();
    out_ptr[0] = my_count;
    for (int i = 0; i < my_count; i++)
        out_ptr[1 + i] = my_scores[i];
    scores_out.release();
}
