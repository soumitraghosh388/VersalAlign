#ifndef TILE_BOTTOM_V3_HH
#define TILE_BOTTOM_V3_HH

#include <adf.h>

// v3: Compute and score collection are DECOUPLED.
// Compute uses batch_in only — no waiting for downstream scores.
// Score output happens after compute is done.

class TileBottomV3 {
private:
    int my_row, my_col;
public:
    TileBottomV3(int row, int col);
    static void registerKernelClass() {
        REGISTER_FUNCTION(TileBottomV3::run);
    }
    void run(
        adf::input_async_buffer<int32_t>& __restrict batch_in,
        adf::output_async_buffer<int32_t>& __restrict scores_out
    );
};

#endif
