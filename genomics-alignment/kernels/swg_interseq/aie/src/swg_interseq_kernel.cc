// Smith-Waterman-Gotoh (SWG) — Inter-sequence SIMD v2 (transposed layout)
//
// swg_interseq_align_pairs(seq1_t, len1, seq2_t, len2, scores):
//   Vectorized: processes BATCH_SIZE=16 pairs simultaneously via SIMD.
//   seq1_t/seq2_t are transposed: index = pos * BATCH_SIZE + pair.
//   Uses aie::load_v<16> at each position to load all 16 pair bases in one
//   instruction, then vectorized max/add for the Gotoh recurrence.
//   scores[BATCH_SIZE] receives the 16 alignment scores on return.
//
//   In aiesim (swg_interseq_kernel): seq1_t holds 16 different query sequences,
//   seq2_t holds 16 different ref sequences.
//   In genomics_server (compute_swg_interseq.hh): seq1_t holds the same query
//   duplicated into all 16 lanes; seq2_t holds 16 different refs.
//
// Single-row in-place DP (same strategy as scalar _v1):
//   M_row[j] / I_row[j] store the previous row's values until overwritten
//   left-to-right.  D is computed purely from the register D_left —
//   no D array is needed (D_above from the previous row is never used in
//   Gotoh's D recurrence, which only depends on the current row's left cell).
//
// Scoring: BLOSUM62 DNA substitution matrix (4×4, A=0 C=1 G=2 T=3)
// Stream protocol: 2-bit encoding (16 bases per int32)

#include <adf.h>
#include "aie_api/aie.hpp"
#include <aie_api/aie_adf.hpp>

#define MAX_SEQ_LEN 128
#define BATCH_SIZE  16
#define GAP_OPEN    -5
#define GAP_EXTEND  -1

// BLOSUM62 DNA 4×4 substitution matrix (A=0 C=1 G=2 T=3)
static const int16_t BLOSUM62_DNA[4][4] = {
    { 4,  0,  0,  0},   // A
    { 0,  9, -3, -1},   // C
    { 0, -3,  6, -2},   // G
    { 0, -1, -2,  5}    // T
};

// Vectorized per-batch SWG — processes BATCH_SIZE=16 pairs in parallel.
//
// Memory layout: single row per matrix (in-place updates).
//   M_row[j*BATCH_SIZE .. (j+1)*BATCH_SIZE - 1] holds M[prev_row, j] until
//   overwritten with M[cur_row, j] as we sweep j left-to-right.
//   Same for I_row.  No D array needed; D is a register-carried value.
//
// On AIE hardware (xchesscc, __chess__ defined): static places M/I rows in
// per-tile data SRAM (~8 KB total, well within 32 KB).
// In x86sim (g++, __chess__ not defined): stack-local allocation avoids the
// shared-static race between concurrent tile instances.
static void swg_interseq_align_pairs(const int8_t* seq1_t, int len1,
                                     const int8_t* seq2_t, int len2,
                                     int16_t scores[BATCH_SIZE]) {
#ifdef __chess__
    alignas(32) static int16_t M_row[(MAX_SEQ_LEN + 1) * BATCH_SIZE];
    alignas(32) static int16_t I_row[(MAX_SEQ_LEN + 1) * BATCH_SIZE];
#else
    alignas(32) int16_t M_row[(MAX_SEQ_LEN + 1) * BATCH_SIZE];
    alignas(32) int16_t I_row[(MAX_SEQ_LEN + 1) * BATCH_SIZE];
#endif

    aie::vector<int16_t, BATCH_SIZE> v_zero       = aie::zeros<int16_t, BATCH_SIZE>();
    aie::vector<int16_t, BATCH_SIZE> v_gap_open   = aie::broadcast<int16_t, BATCH_SIZE>(GAP_OPEN);
    aie::vector<int16_t, BATCH_SIZE> v_gap_extend = aie::broadcast<int16_t, BATCH_SIZE>(GAP_EXTEND);

    // Initialise row 0 of both matrices (boundary condition: all zeros)
    for (int j = 0; j <= len2; j++) {
        aie::store_v(M_row + j * BATCH_SIZE, v_zero);
        aie::store_v(I_row + j * BATCH_SIZE, v_zero);
    }

    aie::vector<int16_t, BATCH_SIZE> max_score = v_zero;

    for (int i = 1; i <= len1; i++) {
        // Load all 16 pairs' base at seq1 position i-1 in one instruction
        aie::vector<int8_t, BATCH_SIZE> base1_i8 = aie::load_v<BATCH_SIZE>(seq1_t + (i - 1) * BATCH_SIZE);

        // Register-carried state for this row sweep (left → right)
        aie::vector<int16_t, BATCH_SIZE> M_diag = v_zero;  // M[i-1, 0] = 0
        aie::vector<int16_t, BATCH_SIZE> M_left = v_zero;  // M[i,   0] = 0
        aie::vector<int16_t, BATCH_SIZE> D_left = v_zero;  // D[i,   0] = 0

        aie::vector<int16_t, BATCH_SIZE> row_max = v_zero;

        for (int j = 1; j <= len2; j++)
            chess_prepare_for_pipelining
            chess_loop_range(1, MAX_SEQ_LEN)
        {
            // Load all 16 pairs' base at seq2 position j-1
            aie::vector<int8_t, BATCH_SIZE> base2_i8 = aie::load_v<BATCH_SIZE>(seq2_t + (j - 1) * BATCH_SIZE);

            // M_row[j] still holds M[i-1, j] (not yet overwritten this sweep)
            aie::vector<int16_t, BATCH_SIZE> M_above = aie::load_v<BATCH_SIZE>(M_row + j * BATCH_SIZE);
            aie::vector<int16_t, BATCH_SIZE> I_above = aie::load_v<BATCH_SIZE>(I_row + j * BATCH_SIZE);

            // Per-lane BLOSUM62 substitution score
            aie::vector<int16_t, BATCH_SIZE> s;
            for (int p = 0; p < BATCH_SIZE; p++) {
                s[p] = BLOSUM62_DNA[(int)base1_i8[p]][(int)base2_i8[p]];
            }

            // Gotoh recurrence
            // I[i,j] = max(M[i-1,j] + gap_open,  I[i-1,j] + gap_extend)
            aie::vector<int16_t, BATCH_SIZE> I_new = aie::max(aie::add(M_above, v_gap_open),
                                                               aie::add(I_above, v_gap_extend));
            // D[i,j] = max(M[i,j-1] + gap_open,  D[i,j-1] + gap_extend)  — pure register
            aie::vector<int16_t, BATCH_SIZE> D_new = aie::max(aie::add(M_left, v_gap_open),
                                                               aie::add(D_left, v_gap_extend));
            // M[i,j] = max(0, M[i-1,j-1]+s, I[i,j], D[i,j])
            aie::vector<int16_t, BATCH_SIZE> M_new = aie::max(aie::max(aie::add(M_diag, s), I_new),
                                                               aie::max(D_new, v_zero));

            // In-place update: overwrite row i-1 values with row i values
            aie::store_v(M_row + j * BATCH_SIZE, M_new);
            aie::store_v(I_row + j * BATCH_SIZE, I_new);

            row_max = aie::max(row_max, M_new);

            // Advance register-carried state
            M_diag = M_above;  // M[i-1,j] becomes M[i-1,j-1] for next j
            D_left = D_new;
            M_left = M_new;
        }

        max_score = aie::max(max_score, row_max);
    }

    for (int p = 0; p < BATCH_SIZE; p++) {
        scores[p] = max_score[p];
    }
}


void swg_interseq_kernel(
    input_stream_int32 *in_seq1,
    input_stream_int32 *in_seq2,
    output_stream_int32 *out_score
) {
    // Transposed sequence buffers: seq[position * BATCH_SIZE + pair]
    // Enables single aie::load_v<16> per position to get all 16 pair bases.
    alignas(32) static int8_t seq1_t[MAX_SEQ_LEN * BATCH_SIZE];
    alignas(32) static int8_t seq2_t[MAX_SEQ_LEN * BATCH_SIZE];

    int32_t batch_count = readincr(in_seq1);

    for (int batch = 0; batch < batch_count; batch++) {
        int32_t lens = readincr(in_seq1);
        int32_t len1 = (lens >> 16) & 0xFFFF;
        int32_t len2 = lens & 0xFFFF;

        // Read 16 seq1 sequences, transpose into seq1_t[pos * BATCH_SIZE + pair]
        for (int p = 0; p < BATCH_SIZE; p++) {
            for (int i = 0; i < (len1 + 15) / 16; i++) {
                int32_t packed = readincr(in_seq1);
                for (int j = 0; j < 16 && (i * 16 + j) < len1; j++) {
                    int pos = i * 16 + j;
                    seq1_t[pos * BATCH_SIZE + p] = (int8_t)((packed >> (j * 2)) & 0x03);
                }
            }
        }

        // Read 16 seq2 sequences, transpose
        int32_t lens2 = readincr(in_seq2);
        for (int p = 0; p < BATCH_SIZE; p++) {
            for (int i = 0; i < (len2 + 15) / 16; i++) {
                int32_t packed = readincr(in_seq2);
                for (int j = 0; j < 16 && (i * 16 + j) < len2; j++) {
                    int pos = i * 16 + j;
                    seq2_t[pos * BATCH_SIZE + p] = (int8_t)((packed >> (j * 2)) & 0x03);
                }
            }
        }

        int16_t scores[BATCH_SIZE];
        swg_interseq_align_pairs(seq1_t, len1, seq2_t, len2, scores);

        for (int p = 0; p < BATCH_SIZE; p += 2) {
            writeincr(out_score, (int32_t)scores[p]);
            writeincr(out_score, (int32_t)scores[p+1]);
        }
    }
}
