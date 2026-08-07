#ifndef COMPUTE_STRIPED_SWG_HH
#define COMPUTE_STRIPED_SWG_HH

// Striped Smith-Waterman compute function for mcast_array.
// Uses build_profile() and sw_striped_align() from sw_striped_kernel_v1.cc.
//
// Batch format (2-bit encoding, 16 bases per int32):
//   batch[0]              = num_refs
//   batch[1]              = query_len
//   batch[2 .. 2+qw-1]   = packed query  (base = (word >> (j*2)) & 0x03, j=0..15)
//   For each ref r:
//     batch[pos]            = ref_len
//     batch[pos+1 .. +rw]   = packed ref  (same packing)
//   where qw = (query_len+15)/16, rw = (ref_len+15)/16
//
// compute_init  — called once per batch; builds the striped query profile.
// compute_score — called per item; aligns ref[item_index] against stored profile.

// Wrap in anonymous namespace so each translation unit gets its own copy,
// avoiding multiple-definition errors in x86sim when both tile_bottom.cc
// and tile_upper.cc include this header.
namespace {
#include "sw_striped_kernel.cc"
}

// Profile passed by reference from the tile's run() stack.
// In x86sim all instances of one kernel class share a single translation
// unit, so a TU-static profile would be overwritten by concurrent tiles
// holding different row data.  On hardware each tile has private SRAM so
// statics are safe, but keeping state here avoids the x86sim race.
struct SwProfile {
    vector<int16, VECLEN> vProfile[AA_LEN * SEG_LEN];  // 4*8*32 B = 1024 B
    int qlen;
};

// Unpack 2-bit encoded bases from packed int32 words into a byte array.
static void ssw_unpack(const int32_t* words, int nwords, int len, int8_t* out) {
    for (int i = 0; i < nwords; i++) {
        int32_t w = words[i];
        for (int j = 0; j < PACKETLEN && (i * PACKETLEN + j) < len; j++)
            out[i * PACKETLEN + j] = (w >> (j * 2)) & 0x03;
    }
}

inline void compute_init(int32_t* batch, int batch_size, SwProfile& prof) {
    int8_t query[MAX_SEQ_LEN];
    int qlen    = batch[1];
    int q_words = (qlen + PACKETLEN - 1) / PACKETLEN;
    ssw_unpack(batch + 2, q_words, qlen, query);

    prof.qlen = qlen;
    build_profile(query, qlen, BLOSUM62_MAT, prof.vProfile);
}

inline int32_t compute_score(int32_t* batch, int batch_size, int item_index,
                              int my_row, int my_col, const SwProfile& prof) {
    int8_t ref[MAX_SEQ_LEN];

    int q_words = (batch[1] + PACKETLEN - 1) / PACKETLEN;
    int pos = 2 + q_words;
    for (int r = 0; r < item_index; r++) {
        int rlen = batch[pos];
        pos += 1 + (rlen + PACKETLEN - 1) / PACKETLEN;
    }

    int rlen    = batch[pos];
    int r_words = (rlen + PACKETLEN - 1) / PACKETLEN;
    ssw_unpack(batch + pos + 1, r_words, rlen, ref);

    return sw_striped_align(ref, rlen, prof.qlen, prof.vProfile);
}

#endif
