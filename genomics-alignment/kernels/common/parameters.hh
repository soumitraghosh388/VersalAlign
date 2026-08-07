#ifndef __PARAMETERS_HH__
#define __PARAMETERS_HH__

// Maximum sequence length supported by DP kernels
// Must match gen_test_data.py and host code
#define MAX_SEQ_LEN 256

// DNA encoding: A=0, C=1, G=2, T=3
#define NUM_BASES 4

// Scoring parameters (affine gap model for SWG)
#define MATCH_SCORE     2
#define MISMATCH_SCORE -1
#define GAP_OPEN       -5
#define GAP_EXTEND     -1

// Linear gap penalty for basic SW
#define GAP_LINEAR     -2

// Debug tracing (0 for performance, 1 for debugging)
// Works in x86sim and aiesim, ignored in HW
#ifndef KERNEL_DEBUG
#define KERNEL_DEBUG 0
#endif

#if KERNEL_DEBUG
#include <stdio.h>
#define TRACE(fmt, ...) printf(fmt "\n", ##__VA_ARGS__)
#else
#define TRACE(fmt, ...)
#endif

#endif // __PARAMETERS_HH__
