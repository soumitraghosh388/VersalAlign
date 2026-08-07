#ifndef COMPUTE_BITAP_HH
#define COMPUTE_BITAP_HH

// Bitap (Wu-Manber shift-OR) compute function for the mcast_array server.
// Uses build_PM() and bitap_align() from bitap_kernel_v1.cc.
//
// Batch format — 2-bit encoding (16 bases per int32), same as compute_striped_swg:
//   batch[0]              = num_refs
//   batch[1]              = query_len  (pattern; must be ≤ BV_WORDS×32 = 128 bp)
//   batch[2 .. 2+qw-1]   = packed query  (base = (word >> (j*2)) & 0x03, j=0..15)
//   For each ref r:
//     batch[pos]            = ref_len
//     batch[pos+1 .. +rw]   = packed ref  (same 2-bit packing)
//   where qw = (query_len+15)/16, rw = (ref_len+15)/16
//
// compute_init  — called once per batch; decodes query/pattern and builds PM[].
// compute_score — called per item; aligns ref[item_index] against stored PM[],
//                 returns minimum edit distance (0..MAX_K) or -1 if > MAX_K.
//
// SwProfile is named to match tile_upper_v3.cc / tile_bottom_v3.cc which
// declare "SwProfile prof" regardless of kernel selection.

// Wrap in anonymous namespace so each translation unit (tile_bottom_v3.cc,
// tile_upper_v3.cc) gets its own copy of bitap_kernel_v1.cc symbols, avoiding
// multiple-definition errors in x86sim — mirrors compute_striped_swg.hh.
namespace {
#include "bitap_kernel.cc"
}

// Profile passed by reference from the tile's run() stack.
// Contains the pattern mask array built by compute_init; used read-only by
// compute_score, so multiple tiles with different patterns don't interfere.
struct SwProfile {
    uint32_t PM[NUM_BASES][BV_WORDS];   // pattern masks (0-bit = match)
    int plen;
    int k;      // min(MAX_K, plen)
};

// Unpack 2-bit encoded bases from packed int32 words into a byte array.
static void bitap_server_unpack(const int32_t* words, int nwords, int len,
                                int8_t* out) {
    for (int i = 0; i < nwords; i++) {
        int32_t w = words[i];
        for (int j = 0; j < 16 && (i * 16 + j) < len; j++)
            out[i * 16 + j] = (w >> (j * 2)) & 0x03;
    }
}

// Decode pattern and build PM[].  Called once per batch before any
// compute_score calls, so PM construction cost is amortized across all refs.
inline void compute_init(int32_t* batch, int batch_size, SwProfile& prof) {
    int8_t pattern[BV_WORDS * 32];     // max 128 bp

    int plen = batch[1];
    if (plen > BV_WORDS * 32) plen = BV_WORDS * 32;    // guard OOB

    int q_words = (plen + 15) / 16;
    bitap_server_unpack(batch + 2, q_words, plen, pattern);

    build_PM(pattern, plen, prof.PM);
    prof.plen = plen;
    prof.k    = (MAX_K < plen) ? MAX_K : plen;
}

// Decode ref at item_index and run Wu-Manber shift-OR against stored PM[].
inline int32_t compute_score(int32_t* batch, int batch_size, int item_index,
                              int my_row, int my_col, const SwProfile& prof) {
    int8_t text[MAX_SEQ_LEN];

    // Walk past query header to reach the ref at item_index
    int q_words = (batch[1] + 15) / 16;
    int pos = 2 + q_words;
    for (int r = 0; r < item_index; r++) {
        int rlen = batch[pos];
        pos += 1 + (rlen + 15) / 16;
    }

    int tlen    = batch[pos];
    int t_words = (tlen + 15) / 16;
    bitap_server_unpack(batch + pos + 1, t_words, tlen, text);

    return bitap_align(text, tlen, prof.plen, prof.PM, prof.k);
}

#endif  // COMPUTE_BITAP_HH
