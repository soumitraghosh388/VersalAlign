#ifndef COMPUTE_SWG_INTERSEQ_HH
#define COMPUTE_SWG_INTERSEQ_HH

// SWG inter-sequence compute adapter for the genomics_server mcast_array.
// Uses swg_interseq_align_pairs() from swg_interseq_kernel.cc, which processes
// BATCH_SIZE=16 refs simultaneously via SIMD.
//
// Key difference from aiesim: here seq1 is ONE query aligned against BATCH_SIZE
// different refs.  To reuse swg_interseq_align_pairs() the query is duplicated
// into all 16 lanes of the transposed seq1_t buffer:
//   seq1_t[pos * BATCH_SIZE + p] = query[pos]  for all p in [0, BATCH_SIZE)
// seq2_t holds BATCH_SIZE different refs:
//   seq2_t[pos * BATCH_SIZE + p] = ref_p[pos]
//
// Interface (unchanged from tile perspective):
//   compute_init()  — decode query; build transposed query_t once per batch
//   compute_score() — on first call (item_index == my_col) pre-computes all
//                     scores for this tile in SIMD batches and caches them;
//                     subsequent calls return cached values
//
// Batch format (2-bit encoding):
//   batch[0]             = num_refs
//   batch[1]             = query_len
//   batch[2 .. 2+qw-1]  = packed query  (base = (word >> (j*2)) & 0x03)
//   For each ref r:
//     batch[pos]           = ref_len
//     batch[pos+1 .. +rw]  = packed ref
//   where qw = (query_len+15)/16, rw = (ref_len+15)/16
//
// SwProfile is named to match tile_upper_v3.cc / tile_bottom_v3.cc.
// In x86sim all instances of one kernel class share a single translation unit,
// so TU-static variables are overwritten by concurrent tiles.  Keeping all
// mutable per-tile state in SwProfile (a local var in each tile's run()) avoids
// the race on both x86sim and hardware.

namespace {
#include "swg_interseq_kernel.cc"
}

struct SwProfile {
    // Query broadcast to all BATCH_SIZE SIMD lanes:
    //   query_t[pos * BATCH_SIZE + p] = query[pos]  for all p
    // alignas(32): aie::load_v<16, int8_t> requires 16-byte alignment; 32 also
    // covers int16_t x16 = 32 bytes used by the DP buffers inside align_pairs().
    alignas(32) int8_t query_t[MAX_SEQ_LEN * BATCH_SIZE];

    // Ref buffer for one SIMD batch: seq2_t[pos*BATCH_SIZE + p] = ref_p[pos].
    // Kept in SwProfile so each tile has its own copy (no x86sim race).
    alignas(32) int8_t seq2_t[MAX_SEQ_LEN * BATCH_SIZE];

    // Pre-computed scores for all refs assigned to this tile.
    // compute_score() fills this on the first call (item_index == my_col) and
    // returns cached values on subsequent calls.
    int32_t scores_cache[64];

    int qlen;
};

inline void compute_init(int32_t* batch, int batch_size, SwProfile& prof) {
    int qlen    = batch[1];
    if (qlen > MAX_SEQ_LEN) qlen = MAX_SEQ_LEN;
    int q_words = (qlen + 15) / 16;
    prof.qlen   = qlen;

    // Decode query and broadcast each base into all BATCH_SIZE lanes.
    // swg_interseq_align_pairs() reads seq1_t[pos*BATCH_SIZE+p] for each pair p,
    // so all lanes must carry the same query base at each position.
    for (int i = 0; i < q_words; i++) {
        int32_t w = batch[2 + i];
        for (int j = 0; j < 16 && (i * 16 + j) < qlen; j++) {
            int8_t base = (int8_t)((w >> (j * 2)) & 0x03);
            int pos = i * 16 + j;
            for (int p = 0; p < BATCH_SIZE; p++)
                prof.query_t[pos * BATCH_SIZE + p] = base;
        }
    }
}

// Internal helper: align up to BATCH_SIZE refs against the query in one SIMD call.
//
// start_item — absolute index of the first ref in this batch (in the full refs list).
//              Subsequent refs are at start_item + num_cols, + 2*num_cols, etc.
// num_cols   — stride between this tile's refs in the full list (= NUM_COLS).
// n          — refs to align: 1 ≤ n ≤ BATCH_SIZE. Unused lanes padded with zeros.
// out_scores — receives n int32_t scores on return.
static void _compute_score_batch(int32_t* batch, int total_items,
                                   int start_item, int num_cols,
                                   SwProfile& prof,
                                   int n, int32_t* out_scores) {
    for (int k = 0; k < MAX_SEQ_LEN * BATCH_SIZE; k++) prof.seq2_t[k] = 0;

    int q_words = (batch[1] + 15) / 16;
    int pos = 2 + q_words;

    for (int r = 0; r < start_item; r++) {
        int rl = batch[pos];
        pos += 1 + (rl + 15) / 16;
    }

    int rlen = 0;
    for (int p = 0; p < n; p++) {
        rlen        = batch[pos];
        int r_words = (rlen + 15) / 16;

        for (int i = 0; i < r_words; i++) {
            int32_t w = batch[pos + 1 + i];
            for (int j = 0; j < 16 && (i * 16 + j) < rlen; j++)
                prof.seq2_t[(i * 16 + j) * BATCH_SIZE + p] = (int8_t)((w >> (j * 2)) & 0x03);
        }
        pos += 1 + r_words;

        if (p < n - 1) {
            for (int skip = 1; skip < num_cols; skip++) {
                int rl = batch[pos];
                pos += 1 + (rl + 15) / 16;
            }
        }
    }

    int16_t scores[BATCH_SIZE];
    swg_interseq_align_pairs(prof.query_t, prof.qlen, prof.seq2_t, rlen, scores);

    for (int p = 0; p < n; p++)
        out_scores[p] = (int32_t)scores[p];
}

// compute_score: called once per ref by tile_bottom_v3.cc / tile_upper_v3.cc.
//
// On the first call for each batch (item_index == my_col), all of this tile's
// refs are pre-computed in SIMD batches of BATCH_SIZE and cached in
// prof.scores_cache[].  Subsequent calls return the cached value for their
// local index.  scores_cache lives in SwProfile so each tile has its own copy.
inline int32_t compute_score(int32_t* batch, int total_items, int item_index,
                              int my_row, int my_col, SwProfile& prof) {
    if (item_index == my_col) {
        int my_count = 0;
        for (int base_i = my_col; base_i < total_items; base_i += NUM_COLS * BATCH_SIZE) {
            int n = 0;
            for (int ii = base_i; ii < total_items && n < BATCH_SIZE; ii += NUM_COLS) n++;
            _compute_score_batch(batch, total_items, base_i, NUM_COLS, prof,
                                 n, &prof.scores_cache[my_count]);
            my_count += n;
        }
    }

    return prof.scores_cache[(item_index - my_col) / NUM_COLS];
}

#endif  // COMPUTE_SWG_INTERSEQ_HH
