// Smith-Waterman-Gotoh (SWG) — Scalar implementation (v3: single-row in-place, no D array)
//
// Restructured from the 2-row ping-pong of earlier versions:
//
//   Old: M_buf[2], I_buf[2], D_buf[2] — six (MAX_SEQ_LEN+1) arrays.
//        D[i,j] = max(M[i,j-1]+go, D[i,j-1]+ge) depends only on the same-row
//        left neighbour, so D_buf stored values that were NEVER READ BACK
//        (D_above was a dead load — Gotoh D recurrence uses D_left register only).
//
//   New: M_row, I_row — two single-row arrays, updated in-place left-to-right.
//        D is a plain register variable (D_left) — no array needed at all.
//
// In-place update (same strategy as swg_interseq_kernel.cc):
//   At position j of row i, M_row[j] still holds M[i-1,j] (not yet overwritten);
//   we read it as M_above, then overwrite with M[i,j].  M_row[0] is initialised
//   to 0 once and never touched by the j=1..len2 loop, so M_diag = 0 always at
//   the start of each row.
//
// Memory: 6 × 129 × 4 = 3,096 bytes → 2 × 129 × 4 = 1,032 bytes (66% reduction)
//
// Scoring: BLOSUM62 DNA substitution matrix (4×4, defined in dna_utils.hh)
// Stream protocol: 2-bit encoding (16 bases per int32)
//
// Recurrence:
//   I[i,j] = max(M[i-1,j] + gap_open,  I[i-1,j] + gap_extend)
//   D[i,j] = max(M[i,j-1] + gap_open,  D[i,j-1] + gap_extend)  — pure register
//   M[i,j] = max(0, M[i-1,j-1]+s, I[i,j], D[i,j])
//   score  = max over all M[i,j]

#include <adf.h>
#include "aie_api/aie.hpp"
#include <aie_api/aie_adf.hpp>
#include "parameters.hh"
#include "dna_utils.hh"

// Core SWG gap-affine alignment for one (seq1, seq2) pair.
// Stack-local DP arrays so this function is safe to call from multiple tiles
// in x86sim (no shared mutable state).  Returns the Smith-Waterman max score.
static int32_t swg_scalar_align(const int8_t* seq1, int len1,
                                 const int8_t* seq2, int len2) {
    int32_t M_row[MAX_SEQ_LEN + 1];
    int32_t I_row[MAX_SEQ_LEN + 1];

    for (int j = 0; j <= len2; j++) {
        M_row[j] = 0;
        I_row[j] = 0;
    }

    int32_t max_score = 0;

    for (int i = 1; i <= len1; i++) {
        int8_t base1 = seq1[i - 1];
        const int16_t* restrict blosum_row = BLOSUM62_DNA[base1];

        int32_t M_diag = 0;
        int32_t M_left = 0;
        int32_t D_left = 0;
        int32_t row_max = 0;

        for (int j = 1; j <= len2; j++)
            chess_prepare_for_pipelining
            chess_loop_range(1, MAX_SEQ_LEN)
        {
            int8_t  base2   = seq2[j - 1];
            int32_t M_above = M_row[j];
            int32_t I_above = I_row[j];
            int32_t s       = (int32_t)blosum_row[base2];

            int32_t i0    = M_above + GAP_OPEN;
            int32_t i1    = I_above + GAP_EXTEND;
            int32_t I_new = (i0 > i1) ? i0 : i1;

            int32_t d0    = M_left + GAP_OPEN;
            int32_t d1    = D_left + GAP_EXTEND;
            int32_t D_new = (d0 > d1) ? d0 : d1;

            aie::vector<int32_t, 4> m_cand;
            m_cand[0] = 0;
            m_cand[1] = M_diag + s;
            m_cand[2] = I_new;
            m_cand[3] = D_new;
            int32_t M_new = aie::reduce_max(m_cand);

            M_row[j] = M_new;
            I_row[j] = I_new;

            if (M_new > row_max) row_max = M_new;

            M_diag = M_above;
            D_left = D_new;
            M_left = M_new;
        }

        max_score = (row_max > max_score) ? row_max : max_score;
    }

    return max_score;
}

void swg_kernel(
    input_stream_int32 *in_seq1,
    input_stream_int32 *in_seq2,
    output_stream_int32 *out_score
) {
    static int8_t seq1[MAX_SEQ_LEN];
    static int8_t seq2[MAX_SEQ_LEN];

    int32_t pair_count = readincr(in_seq1);
    TRACE("SWG: pair_count=%d", pair_count);

    for (int pair = 0; pair < pair_count; pair++) {
        int32_t lens = readincr(in_seq1);
        int32_t len1 = (lens >> 16) & 0xFFFF;
        int32_t len2 = lens & 0xFFFF;
        TRACE("SWG: pair=%d len1=%d len2=%d", pair, len1, len2);

        for (int i = 0; i < (len1 + 15) / 16; i++) {
            int32_t packed = readincr(in_seq1);
            int remaining = ((i * 16 + 16) <= len1) ? 16 : (len1 - i * 16);
            unpack_bases_2bit(packed, &seq1[i * 16], remaining);
        }

        int32_t lens2 = readincr(in_seq2);
        for (int i = 0; i < (len2 + 15) / 16; i++) {
            int32_t packed = readincr(in_seq2);
            int remaining = ((i * 16 + 16) <= len2) ? 16 : (len2 - i * 16);
            unpack_bases_2bit(packed, &seq2[i * 16], remaining);
        }

        int32_t max_score = swg_scalar_align(seq1, len1, seq2, len2);
        TRACE("SWG: pair=%d score=%d", pair, max_score);
        writeincr(out_score, max_score);
    }
}
