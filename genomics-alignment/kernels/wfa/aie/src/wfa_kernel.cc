// Wavefront Alignment (WFA) Kernel for AMD Versal AIE
//
// Implements the WFA algorithm (Marco-Sola et al., Bioinformatics 2021).
// Expands wavefronts score-by-score until the target cell is reached.
// Complexity: O(s^2 + s*n) where s = alignment penalty.
// For similar sequences (low s), much faster than O(mn) DP.
//
// Substitution penalty: fixed WFA_MISMATCH = 4 for any mismatched base pair.
// Exact matches cost 0 and are handled by the extend step.
//
// Gap-affine penalties: gap_open=5, gap_extend=1  (matches wfa_align in generate_test_data.py)
// Returns the alignment penalty (lower = better; -1 if > MAX_SCORE).
//
// Stream protocol: 2-bit encoding (16 bases per int32).
//   seq1 stream: [pair_count], per pair: [len1<<16|len2], 2-bit packed bases...
//   seq2 stream: per pair: [sync_word], 2-bit packed bases...
//   output stream: one alignment penalty per pair

#include <adf.h>
#include "aie_api/aie.hpp"
#include <aie_api/aie_adf.hpp>
#include "dna_utils.hh"

// Gap-affine penalties
#define WFA_GAP_OPEN  5
#define WFA_GAP_EXT   1

// Fixed mismatch penalty (any mismatched base pair)
#define WFA_MISMATCH 4

// Maximum alignment penalty to search
#define MAX_SCORE 60// (MAX_LEN+MAX_LEN)*WFA_MISMATCH=1024

// Wavefront array width: diagonals k = -MAX_SCORE .. +MAX_SCORE
#define MAX_WF_WIDTH (2 * MAX_SCORE + 1)

// Sentinel: no reachable wavefront on this diagonal
#define WF_NULL (-999999)

// Circular buffer depth.
// Gap-open lookback is WFA_GAP_OPEN+WFA_GAP_EXT (matches WFA2-lib's
// wavefront_compute_fetch_input(): gap_open1 = score - gap_opening1 - gap_extension1,
// i.e. a gap of length L costs open + extend*L, not open + extend*(L-1)).
// Max lookback = max(WFA_GAP_OPEN+WFA_GAP_EXT, WFA_MISMATCH) = max(6, 4) = 6  → 7 slots.
#define WF_LOOKBACK 7

// Run WFA gap-affine alignment for one (seq1, seq2) pair.
// All wavefront arrays are stack-local so this function is safe to call from
// multiple tiles in x86sim (no shared mutable state).
// Returns alignment penalty 0..MAX_SCORE, or -1 if penalty exceeds MAX_SCORE.
static int32_t wfa_align(const int8_t* seq1, int plen,
                          const int8_t* seq2, int tlen) {
    // Stack-local wavefront circular buffer (~4.4 KB)
    int16_t wf_m[WF_LOOKBACK][MAX_WF_WIDTH];
    int16_t wf_i[WF_LOOKBACK][MAX_WF_WIDTH];
    int16_t wf_d[WF_LOOKBACK][MAX_WF_WIDTH];
    int16_t slot_score[WF_LOOKBACK];

    int target_k   = plen - tlen;
    int target_off = plen;
    int center     = MAX_SCORE;

    for (int slot = 0; slot < WF_LOOKBACK; slot++) {
        slot_score[slot] = -1;
        for (int d = 0; d < MAX_WF_WIDTH; d++) {
            wf_m[slot][d] = WF_NULL;
            wf_i[slot][d] = WF_NULL;
            wf_d[slot][d] = WF_NULL;
        }
    }

    slot_score[0] = 0;
    wf_m[0][center] = 0;

    {
        int16_t off = wf_m[0][center];
        int i = off, j = off;
        while (i < plen && j < tlen && seq1[i] == seq2[j]) { off++; i++; j++; }
        wf_m[0][center] = off;
    }

    int32_t result = -1;
    if (wf_m[0][center + target_k] >= target_off) result = 0;

    for (int s = 1; s <= MAX_SCORE && result < 0; s++) {
        int cur_slot = s % WF_LOOKBACK;

        for (int d = 0; d < MAX_WF_WIDTH; d++) {
            wf_m[cur_slot][d] = WF_NULL;
            wf_i[cur_slot][d] = WF_NULL;
            wf_d[cur_slot][d] = WF_NULL;
        }
        slot_score[cur_slot] = s;

        int s_go = s - WFA_GAP_OPEN - WFA_GAP_EXT;
        int s_ge = s - WFA_GAP_EXT;
        int slot_go = (s_go >= 0) ? (s_go % WF_LOOKBACK) : -1;
        int slot_ge = (s_ge >= 0) ? (s_ge % WF_LOOKBACK) : -1;
        if (slot_go >= 0 && slot_score[slot_go] != s_go) slot_go = -1;
        if (slot_ge >= 0 && slot_score[slot_ge] != s_ge) slot_ge = -1;

        int lo = center - s; if (lo < 0) lo = 0;
        int hi = center + s; if (hi >= MAX_WF_WIDTH) hi = MAX_WF_WIDTH - 1;

        // ── Insertion (gap in seq1) ──────────────────────────────────────────
        for (int d = lo; d <= hi; d++) {
            int16_t val = WF_NULL;
            if (slot_go >= 0 && d + 1 < MAX_WF_WIDTH) {
                int16_t v = wf_m[slot_go][d + 1];
                if (v != WF_NULL) { v += 1; if (v > val) val = v; }
            }
            if (slot_ge >= 0 && d + 1 < MAX_WF_WIDTH) {
                int16_t v = wf_i[slot_ge][d + 1];
                if (v != WF_NULL) { v += 1; if (v > val) val = v; }
            }
            wf_i[cur_slot][d] = val;
        }

        // ── Deletion (gap in seq2) ───────────────────────────────────────────
        for (int d = lo; d <= hi; d++) {
            int16_t val = WF_NULL;
            if (slot_go >= 0 && d - 1 >= 0) {
                int16_t v = wf_m[slot_go][d - 1];
                if (v != WF_NULL && v > val) val = v;
            }
            if (slot_ge >= 0 && d - 1 >= 0) {
                int16_t v = wf_d[slot_ge][d - 1];
                if (v != WF_NULL && v > val) val = v;
            }
            wf_d[cur_slot][d] = val;
        }

        // ── Match wavefront ──────────────────────────────────────────────────
        for (int d = lo; d <= hi; d++) {
            int16_t val = WF_NULL;
            int k = d - center;

            {
                int s_sub = s - WFA_MISMATCH;
                if (s_sub >= 0) {
                    int slot_sub = s_sub % WF_LOOKBACK;
                    if (slot_score[slot_sub] == s_sub) {
                        int16_t v = wf_m[slot_sub][d];
                        if (v != WF_NULL) {
                            int ip = (int)v, jp = ip - k;
                            if (ip >= 0 && ip < plen && jp >= 0 && jp < tlen
                                    && seq1[ip] != seq2[jp]) {
                                int16_t nv = v + 1;
                                if (nv > val) val = nv;
                            }
                        }
                    }
                }
            }
            { int16_t v = wf_i[cur_slot][d]; if (v != WF_NULL && v > val) val = v; }
            { int16_t v = wf_d[cur_slot][d]; if (v != WF_NULL && v > val) val = v; }
            wf_m[cur_slot][d] = val;
        }

        // ── Extend along exact matches ───────────────────────────────────────
        for (int d = lo; d <= hi; d++) {
            int16_t off = wf_m[cur_slot][d];
            if (off == WF_NULL) continue;
            int k = d - center, i = off, j = off - k;
            while (i < plen && j < tlen && i >= 0 && j >= 0 && seq1[i] == seq2[j]) {
                off++; i++; j++;
            }
            wf_m[cur_slot][d] = off;
        }

        // ── Check termination ────────────────────────────────────────────────
        int target_d = center + target_k;
        if (target_d >= lo && target_d <= hi && wf_m[cur_slot][target_d] >= target_off)
            result = s;
    }
    return result;
}

void wfa_kernel(
    input_stream_int32 *in_seq1,
    input_stream_int32 *in_seq2,
    output_stream_int32 *out_score
) {
    static int8_t seq1[MAX_SEQ_LEN];
    static int8_t seq2[MAX_SEQ_LEN];

    int32_t pair_count = readincr(in_seq1);

    for (int pair = 0; pair < pair_count; pair++) {
        int32_t lens = readincr(in_seq1);
        int32_t plen = (lens >> 16) & 0xFFFF;
        int32_t tlen = lens & 0xFFFF;

        int pat_words = (plen + 15) / 16;
        for (int i = 0; i < pat_words; i++) {
            int32_t packed = readincr(in_seq1);
            for (int j = 0; j < 16 && (i * 16 + j) < plen; j++)
                seq1[i * 16 + j] = (packed >> (j * 2)) & 0x03;
        }

        int32_t lens2 = readincr(in_seq2);
        int txt_words = (tlen + 15) / 16;
        for (int i = 0; i < txt_words; i++) {
            int32_t packed = readincr(in_seq2);
            for (int j = 0; j < 16 && (i * 16 + j) < tlen; j++)
                seq2[i * 16 + j] = (packed >> (j * 2)) & 0x03;
        }

        writeincr(out_score, wfa_align(seq1, plen, seq2, tlen));
    }
}
