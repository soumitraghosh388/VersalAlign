#ifndef COMPUTE_DUMMY_HH
#define COMPUTE_DUMMY_HH

// Dummy compute: returns value + tile-specific offset.
// For scaffolding and data path verification.

#include <cstdint>

inline void compute_init(int32_t* /*batch*/, int /*batch_size*/) {
    // No initialization needed for dummy
}

inline int32_t compute_score(int32_t* batch, int batch_size, int item_index,
                              int my_row, int my_col) {
    int32_t val = batch[1 + item_index];  // skip header word
    return val + (my_col + 1) * 100 + (my_row + 1);
}

#endif
