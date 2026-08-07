#ifndef TILE_UPPER_V3_HH
#define TILE_UPPER_V3_HH

#include <adf.h>

// v3: Compute FIRST using only batch_in (no waiting for downstream).
// THEN collect downstream scores and merge with own scores.
// This allows all rows to compute in PARALLEL.

class TileUpperV3 {
private:
    int my_row, my_col;
public:
    TileUpperV3(int row, int col);
    static void registerKernelClass() {
        REGISTER_FUNCTION(TileUpperV3::run);
    }
    void run(
        adf::input_async_buffer<int32_t>& __restrict batch_in,
        adf::input_async_buffer<int32_t>& __restrict scores_in,
        adf::output_async_buffer<int32_t>& __restrict scores_out
    );
};

#endif
