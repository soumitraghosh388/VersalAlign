// Bitap (Shift-OR) Edit Distance Kernel for AMD Versal AIE
//
// Implements Wu-Manber's approximate string matching for edit distance.
// Uses k+1 bit-vectors with only shifts, ANDs, and ORs.
// Multi-word shifts have NO carry chain issues (unlike BPM's addition).
//
// This version: MAX_K configurable, fixed 4-word (128-bit) bit-vectors for up to 128bp.
// Levels 0-5 are fully unrolled; levels 6..MAX_K use a scalar loop.
//
// Note: BLOSUM62 does not apply here — pattern masks are binary (0=match, 1=mismatch).
// All substitutions cost 1 edit regardless of base identity (Wu-Manber model).
//
// Stream protocol: 2-bit encoding (16 bases per int32), matches swg_scalar.
//
// Complexity: O(k * n) where k = max errors, n = text length.
// For Illumina (100bp, k=10): ~72 scalar ops per text char, ~7200 total.

#include <adf.h>
#include "aie_api/aie.hpp"
#include <aie_api/aie_adf.hpp>
#include "dna_utils.hh"

#define MAX_K 5
#define BV_WORDS 4  // 4 x 32 = 128 bits max pattern length

// Inlined 4-word left shift by 1 bit
#define BV_SHL1(dst, src) do { \
    (dst)[3] = ((src)[3] << 1) | ((src)[2] >> 31); \
    (dst)[2] = ((src)[2] << 1) | ((src)[1] >> 31); \
    (dst)[1] = ((src)[1] << 1) | ((src)[0] >> 31); \
    (dst)[0] = (src)[0] << 1; \
} while(0)

// Inlined 4-word OR: dst = a | b
#define BV_OR(dst, a, b) do { \
    (dst)[0] = (a)[0] | (b)[0]; \
    (dst)[1] = (a)[1] | (b)[1]; \
    (dst)[2] = (a)[2] | (b)[2]; \
    (dst)[3] = (a)[3] | (b)[3]; \
} while(0)

// Inlined 4-word AND: dst = a & b
#define BV_AND(dst, a, b) do { \
    (dst)[0] = (a)[0] & (b)[0]; \
    (dst)[1] = (a)[1] & (b)[1]; \
    (dst)[2] = (a)[2] & (b)[2]; \
    (dst)[3] = (a)[3] & (b)[3]; \
} while(0)

// Inlined 4-word copy
#define BV_COPY(dst, src) do { \
    (dst)[0] = (src)[0]; \
    (dst)[1] = (src)[1]; \
    (dst)[2] = (src)[2]; \
    (dst)[3] = (src)[3]; \
} while(0)

// Inlined: dst = (src << 1) | pm
#define BV_SHL1_OR(dst, src, pm) do { \
    (dst)[3] = ((src)[3] << 1) | ((src)[2] >> 31) | (pm)[3]; \
    (dst)[2] = ((src)[2] << 1) | ((src)[1] >> 31) | (pm)[2]; \
    (dst)[1] = ((src)[1] << 1) | ((src)[0] >> 31) | (pm)[1]; \
    (dst)[0] = ((src)[0] << 1) | (pm)[0]; \
} while(0)

// Process one error level: R[d] = ((R[d]<<1)|pm) & (old[d-1]<<1) & old[d-1] & (new[d-1]<<1)
// All temps computed inline
#define BITAP_LEVEL(Rd, old_Rd, pm, old_Rd1, new_Rd1, t1, t2, t3) do { \
    BV_SHL1_OR(t1, old_Rd, pm);   /* t1 = (R[d]<<1) | PM  (match extension) */ \
    BV_SHL1(t2, old_Rd1);          /* t2 = old_R[d-1] << 1 (substitution)    */ \
    BV_AND(t1, t1, t2);            /* t1 &= substitution                     */ \
    BV_AND(t1, t1, old_Rd1);       /* t1 &= old_R[d-1]     (insertion)       */ \
    BV_SHL1(t3, new_Rd1);          /* t3 = R[d-1] << 1      (deletion)       */ \
    BV_AND(Rd, t1, t3);            /* R[d] = t1 & deletion                   */ \
} while(0)

// Build pattern mask array from a decoded pattern byte sequence.
// PM[base][word]: bit clear at position pos means pattern[pos] == base.
// Equivalent to the inline init block inside bitap_kernel(), extracted for
// reuse by compute_bitap.hh (analogous to build_profile() in sw_striped).
static void build_PM(const int8_t* pattern, int plen,
                     uint32_t PM[NUM_BASES][BV_WORDS]) {
    for (int c = 0; c < NUM_BASES; c++)
        PM[c][0] = PM[c][1] = PM[c][2] = PM[c][3] = 0xFFFFFFFF;
    for (int pos = 0; pos < plen; pos++)
        PM[(int)pattern[pos]][pos >> 5] &= ~(1U << (pos & 31));
}

// Run Wu-Manber shift-OR for one (text, PM) pair.
// R/old_R/t1/t2/t3 are stack-local so this function is safe to call from
// multiple tiles in x86sim (no shared mutable state).
// Returns minimum edit distance 0..k, or -1 if no alignment within k errors.
static int32_t bitap_align(const int8_t* text, int tlen, int plen,
                            const uint32_t PM[NUM_BASES][BV_WORDS], int k) {
    uint32_t R    [MAX_K + 1][BV_WORDS];
    uint32_t old_R[MAX_K + 1][BV_WORDS];
    uint32_t t1[BV_WORDS], t2[BV_WORDS], t3[BV_WORDS];

    for (int d = 0; d <= k; d++) {
        R[d][0] = R[d][1] = R[d][2] = R[d][3] = 0xFFFFFFFF;
        for (int i = 0; i < d && i < plen; i++)
            R[d][i >> 5] &= ~(1U << (i & 31));
    }

    for (int j = 0; j < tlen; j++) {
        const uint32_t* pm = PM[(int)text[j]];
        for (int d = 0; d <= k; d++) BV_COPY(old_R[d], R[d]);

        BV_SHL1_OR(R[0], old_R[0], pm);
        if (k >= 1)  BITAP_LEVEL(R[1],  old_R[1],  pm, old_R[0],  R[0],  t1, t2, t3);
        if (k >= 2)  BITAP_LEVEL(R[2],  old_R[2],  pm, old_R[1],  R[1],  t1, t2, t3);
        if (k >= 3)  BITAP_LEVEL(R[3],  old_R[3],  pm, old_R[2],  R[2],  t1, t2, t3);
        if (k >= 4)  BITAP_LEVEL(R[4],  old_R[4],  pm, old_R[3],  R[3],  t1, t2, t3);
        if (k >= 5)  BITAP_LEVEL(R[5],  old_R[5],  pm, old_R[4],  R[4],  t1, t2, t3);
        if (k >= 6)  BITAP_LEVEL(R[6],  old_R[6],  pm, old_R[5],  R[5],  t1, t2, t3);
        if (k >= 7)  BITAP_LEVEL(R[7],  old_R[7],  pm, old_R[6],  R[6],  t1, t2, t3);
        if (k >= 8)  BITAP_LEVEL(R[8],  old_R[8],  pm, old_R[7],  R[7],  t1, t2, t3);
        if (k >= 9)  BITAP_LEVEL(R[9],  old_R[9],  pm, old_R[8],  R[8],  t1, t2, t3);
        if (k >= 10) BITAP_LEVEL(R[10], old_R[10], pm, old_R[9],  R[9],  t1, t2, t3);
        for (int d = 11; d <= k; d++) {
            BITAP_LEVEL(R[d], old_R[d], pm, old_R[d-1], R[d-1], t1, t2, t3);
        }
    }

    int check_word     = (plen - 1) >> 5;
    uint32_t check_bit = 1U << ((plen - 1) & 31);
    for (int d = 0; d <= k; d++) {
        if ((R[d][check_word] & check_bit) == 0) return d;
    }
    return -1;
}

void bitap_kernel(
    input_stream_int32 *in_seq1,
    input_stream_int32 *in_seq2,
    output_stream_int32 *out_score
) {
    static int8_t    pattern[MAX_SEQ_LEN];
    static int8_t    text[MAX_SEQ_LEN];
    static uint32_t  PM[NUM_BASES][BV_WORDS];

    int32_t pair_count = readincr(in_seq1);

    for (int pair = 0; pair < pair_count; pair++) {
        int32_t lens = readincr(in_seq1);
        int32_t plen = (lens >> 16) & 0xFFFF;
        int32_t tlen = lens & 0xFFFF;

        // Read pattern — 2-bit encoding: 16 bases per int32
        int pattern_words = (plen + 15) / 16;
        for (int i = 0; i < pattern_words; i++) {
            int32_t packed = readincr(in_seq1);
            for (int j = 0; j < 16 && (i * 16 + j) < plen; j++)
                pattern[i * 16 + j] = (packed >> (j * 2)) & 0x03;
        }

        build_PM(pattern, plen, PM);

        // Read text — 2-bit encoding: 16 bases per int32
        int32_t lens2 = readincr(in_seq2);
        int text_words = (tlen + 15) / 16;
        for (int i = 0; i < text_words; i++) {
            int32_t packed = readincr(in_seq2);
            for (int j = 0; j < 16 && (i * 16 + j) < tlen; j++)
                text[i * 16 + j] = (packed >> (j * 2)) & 0x03;
        }

        int k = MAX_K;
        if (k > plen) k = plen;

        writeincr(out_score, bitap_align(text, tlen, plen, PM, k));
    }
}
