#ifndef COMPUTE_WFA_HH
#define COMPUTE_WFA_HH

// WFA gap-affine compute function for the mcast_array server.
// Uses wfa_align() from wfa_kernel_v1.cc.
//
// Batch format (2-bit encoding, 16 bases per int32, same as compute_striped_swg):
//   batch[0]             = num_refs
//   batch[1]             = query_len
//   batch[2 .. 2+qw-1]  = packed query  (base = (word >> (j*2)) & 0x03, j=0..15)
//   For each ref r:
//     batch[pos]           = ref_len
//     batch[pos+1 .. +rw]  = packed ref  (same 2-bit packing)
//   where qw = (query_len+15)/16, rw = (ref_len+15)/16
//
// compute_init  — called once per batch; decodes query into SwProfile.
// compute_score — called per item; decodes ref[item_index] and calls wfa_align().
//
// SwProfile is named to match tile_upper_v3.cc / tile_bottom_v3.cc which
// declare "SwProfile prof" regardless of kernel selection.

// Wrap in anonymous namespace so each translation unit (tile_bottom_v3.cc,
// tile_upper_v3.cc) gets its own copy of wfa_kernel_v1.cc symbols, avoiding
// multiple-definition errors in x86sim — mirrors compute_striped_swg.hh.
namespace {
#include "wfa_kernel.cc"
}

// Profile holds the decoded query for the current batch.
// wfa_align() needs both sequences at call time, so only the query is cached.
struct SwProfile {
    int8_t query[MAX_SEQ_LEN];
    int qlen;
};

static void wfa_server_unpack(const int32_t* words, int nwords, int len, int8_t* out) {
    for (int i = 0; i < nwords; i++) {
        int32_t w = words[i];
        for (int j = 0; j < 16 && (i * 16 + j) < len; j++)
            out[i * 16 + j] = (w >> (j * 2)) & 0x03;
    }
}

inline void compute_init(int32_t* batch, int batch_size, SwProfile& prof) {
    int qlen    = batch[1];
    if (qlen > MAX_SEQ_LEN) qlen = MAX_SEQ_LEN;
    int q_words = (qlen + 15) / 16;
    wfa_server_unpack(batch + 2, q_words, qlen, prof.query);
    prof.qlen   = qlen;
}

inline int32_t compute_score(int32_t* batch, int batch_size, int item_index,
                              int my_row, int my_col, const SwProfile& prof) {
    int8_t ref[MAX_SEQ_LEN];

    int q_words = (batch[1] + 15) / 16;
    int pos = 2 + q_words;
    for (int r = 0; r < item_index; r++) {
        int rlen = batch[pos];
        pos += 1 + (rlen + 15) / 16;
    }

    int rlen    = batch[pos];
    int r_words = (rlen + 15) / 16;
    wfa_server_unpack(batch + pos + 1, r_words, rlen, ref);

    return wfa_align(prof.query, prof.qlen, ref, rlen);
}

#endif  // COMPUTE_WFA_HH
