// v3 Upper tile: compute FIRST, collect scores AFTER.
//
// Phase 1: Acquire batch_in → compute all my scores → release batch_in
//          (scores_in NOT acquired yet — so this tile doesn't wait for below)
//
// Phase 2: Acquire scores_in → read downstream scores
//          Acquire scores_out → write downstream + my scores → release both
//
// All rows compute Phase 1 in PARALLEL (no inter-row dependency).
// Phase 2 is sequential bottom-to-top but very fast (just copying scores).

#include "tile_upper_v3.hh"
#include <aie_api/aie.hpp>
#include "kernel_config.h"
#include "compute/compute.hh"

TileUpperV3::TileUpperV3(int row, int col) : my_row(row), my_col(col) {}

void TileUpperV3::run(
    adf::input_async_buffer<int32_t>& __restrict batch_in,
    adf::input_async_buffer<int32_t>& __restrict scores_in,
    adf::output_async_buffer<int32_t>& __restrict scores_out
) {
    // Phase 1: COMPUTE (only batch_in — all rows run in parallel)
    batch_in.acquire();
    auto in_ptr = (int32_t*)batch_in.data();
    int32_t total_items = in_ptr[0];

    SwProfile prof;
    compute_init(in_ptr, total_items, prof);

    int32_t my_scores[64];  // max refs per tile
    int my_count = 0;
    for (int i = my_col; i < total_items; i += NUM_COLS) {
        my_scores[my_count] = compute_score(in_ptr, total_items, i, my_row, my_col, prof);
        my_count++;
    }
    batch_in.release();

    // Phase 2: COLLECT SCORES (sequential bottom-to-top, but fast)
    scores_in.acquire();
    auto sin_ptr = (int32_t*)scores_in.data();
    int32_t downstream_count = sin_ptr[0];

    scores_out.acquire();
    auto out_ptr = (int32_t*)scores_out.data();

    // Downstream scores first
    for (int i = 0; i < downstream_count; i++)
        out_ptr[1 + i] = sin_ptr[1 + i];

    // Append my scores
    for (int i = 0; i < my_count; i++)
        out_ptr[1 + downstream_count + i] = my_scores[i];

    out_ptr[0] = downstream_count + my_count;

    scores_in.release();
    scores_out.release();
}
