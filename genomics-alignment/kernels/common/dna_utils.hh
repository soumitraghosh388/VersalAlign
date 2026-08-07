#ifndef __DNA_UTILS_HH__
#define __DNA_UTILS_HH__

#include "parameters.hh"

// ---- BLOSUM62 for DNA (4×4, A=0 C=1 G=2 T=3) --------------------------------
//
// Scores based on nucleotide substitution rates:
//   Match      : +2  (same base)
//   Transition : -1  (A↔G purine swap, C↔T pyrimidine swap — more frequent)
//   Transversion: -2 (purine↔pyrimidine — less frequent, larger penalty)
//
//   Transitions (rate ~2× higher than transversions in human genome):
//     A↔G (row 0 col 2, row 2 col 0)
//     C↔T (row 1 col 3, row 3 col 1)
//
// Replaces flat MATCH_SCORE/MISMATCH_SCORE in the SWG inner loop.

static const int16_t BLOSUM62_DNA[4][4] = {
//   A   C   G   T
    { 4, 0, 0, 0},   // A: match=4, transv C=0, transit G=0, transv T=0
    {0,  9, -3, -1},   // C: transv A=0, match=9, transv G=-3, transit T=-1
    {0, -3,  6, -2},   // G: transit A=0, transv C=-3, match=6, transv T=-2
    {0, -1, -2,  5}    // T: transv A=0, transit C=-1, transv G=-2, match=5
};

// Inline lookup — use cached row pointer in loops for pipelining efficiency
inline int16_t score_blosum62(int8_t base1, int8_t base2) {
    return BLOSUM62_DNA[base1][base2];
}

// Legacy flat scoring (kept for other kernels that still use it)
inline int16_t score_match(int8_t base1, int8_t base2) {
    return (base1 == base2) ? MATCH_SCORE : MISMATCH_SCORE;
}

// ---- Max helpers ------------------------------------------------------------
inline int16_t max2(int16_t a, int16_t b) {
    return (a > b) ? a : b;
}

inline int16_t max3(int16_t a, int16_t b, int16_t c) {
    int16_t m = (a > b) ? a : b;
    return (m > c) ? m : c;
}

inline int16_t max4(int16_t a, int16_t b, int16_t c, int16_t d) {
    int16_t m1 = (a > b) ? a : b;
    int16_t m2 = (c > d) ? c : d;
    return (m1 > m2) ? m1 : m2;
}

// ---- Base unpacking ---------------------------------------------------------
// 2-bit encoding: 16 bases per int32 (A=0,C=1,G=2,T=3)
inline void unpack_bases_2bit(int32_t packed, int8_t* bases, int count) {
    for (int j = 0; j < count; j++)
        bases[j] = (packed >> (j * 2)) & 0x03;
}

// 8-bit encoding: 4 bases per int32 (legacy, kept for other kernels)
inline void unpack_bases(int32_t packed, int8_t* bases, int count) {
    for (int j = 0; j < count; j++)
        bases[j] = (packed >> (j * 8)) & 0x03;
}

#endif // __DNA_UTILS_HH__
