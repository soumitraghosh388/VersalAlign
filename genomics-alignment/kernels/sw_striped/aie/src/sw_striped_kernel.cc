// Smith-Waterman — Striped SIMD (Farrar 2007)
//
// Ported from ~/versal-acap/genome_sequence/SWG/aie/src_stripped_st_sk_10ref/
// Key patterns from the working reference:
//   - vH.push(0) for vector shift (native AIE API)
//   - Unsigned subtraction model: saturating_subu(H, gap) = max(H - gap, 0)
//   - Gap penalties are positive values, subtracted (not negative values added)
//   - vector<int16, 16> stored directly in arrays
//
// BLOSUM62 DNA scoring, 2-bit base encoding (16 bases per int32).
// VECLEN=16, max query=128, max ref=128

#include <adf.h>
#include "aie_api/aie.hpp"
#include <aie_api/aie_adf.hpp>

using namespace aie;

#define VECLEN    16
#define MAX_SEQ_LEN 128
#define SEG_LEN ((MAX_SEQ_LEN + VECLEN - 1) / VECLEN)  // 8 segments max
#define AA_LEN    4
#define PACKETLEN 16   // 2-bit encoding: 16 bases per int32

#define GAP_OPEN 5    // positive, will be subtracted (matches BLOSUM62 convention)
#define GAP_EXT  1    // positive, will be subtracted

// Out-of-bounds profile fill (pad segments beyond readLen)
const int def_mismatch = 10;

// BLOSUM62 DNA 4×4 substitution matrix (A=0 C=1 G=2 T=3), row-major
static const int8_t BLOSUM62_MAT[AA_LEN * AA_LEN] = {
     4,  0,  0,  0,   // A: match=4, transv C=0, transit G=0, transv T=0
     0,  9, -3, -1,   // C: transv A=0, match=9, transv G=-3, transit T=-1
     0, -3,  6, -2,   // G: transit A=0, transv C=-3, match=6, transv T=-2
     0, -1, -2,  5    // T: transv A=0, transit C=-1, transv G=-2, match=5
};

#define saturating_subu(a, b) max(vector<int16, VECLEN>(a - b), zeros<int16, VECLEN>())

// Build query profile in striped order.
// vProfile[nt * segLen + seg] holds the BLOSUM62 score for nucleotide nt, segment seg.
void build_profile(const int8_t* read, int readLen,
                   const int8_t* mat,
                   vector<int16, VECLEN>* vProfile) {
    int segLen = (readLen + VECLEN - 1) / VECLEN;

    for (int nt = 0; nt < AA_LEN; nt++) {
        for (int i = 0; i < segLen; i++) {
            vector<int16, VECLEN> vals = zeros<int16, VECLEN>();
            int j = i;
            for (int s = 0; s < VECLEN; s++) {
                vals[s] = (j >= readLen) ? -def_mismatch : mat[nt * AA_LEN + read[j]];
                j += segLen;
            }
            vProfile[nt * segLen + i] = vals;
        }
    }
}

// Striped SW alignment (single pair)
int32_t sw_striped_align(const int8_t* ref, int refLen,
                         int readLen,
                         const vector<int16, VECLEN>* vProfile) {
    int segLen = (readLen + VECLEN - 1) / VECLEN;

    vector<int16, VECLEN> vZero = zeros<int16, VECLEN>();
    vector<int16, VECLEN> vGapO = broadcast<int16, VECLEN>((int16)GAP_OPEN);
    vector<int16, VECLEN> vGapE = broadcast<int16, VECLEN>((int16)GAP_EXT);

    vector<int16, VECLEN> pvHStoreA[SEG_LEN];
    vector<int16, VECLEN> pvHLoadA[SEG_LEN];
    vector<int16, VECLEN> pvEA[SEG_LEN];

    for (int i = 0; i < segLen; i++) {
        pvHStoreA[i] = vZero;
        pvHLoadA[i] = vZero;
        pvEA[i] = vZero;
    }

    vector<int16, VECLEN>* pvHStore = pvHStoreA;
    vector<int16, VECLEN>* pvHLoad = pvHLoadA;
    vector<int16, VECLEN>* pvE = pvEA;

    vector<int16, VECLEN> vMaxScore = vZero;

    for (int i = 0; i < refLen; i++) {
        vector<int16, VECLEN> e, vF = vZero;
        vector<int16, VECLEN> vH = pvHStore[segLen - 1];
        vH.push(0);  // Shift right, insert 0 at position 0

        // Swap H buffers
        vector<int16, VECLEN>* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        vector<int16, VECLEN> vMaxColumn = vZero;
        const vector<int16, VECLEN>* vP = vProfile + ref[i] * segLen;

        // Process each segment
        for (int j = 0; j < segLen; j++) {
            // Add profile score (saturating)
            vH = min(max(vector<int16, VECLEN>(vH + vP[j]),
                         broadcast<int16, VECLEN>((int16)-32768)),
                     broadcast<int16, VECLEN>((int16)32767));

            // H = max(H, E, F)
            e = pvE[j];
            vH = max(vH, e);
            vH = max(vH, vF);
            vMaxColumn = max(vMaxColumn, vH);

            // Store H
            pvHStore[j] = vH;

            // Update E: E = max(H - gapO, E - gapE)
            vH = saturating_subu(vH, vGapO);
            e = saturating_subu(e, vGapE);
            e = max(e, vH);
            pvE[j] = e;

            // Update F: F = max(H - gapO, F - gapE)
            vF = saturating_subu(vF, vGapE);
            vF = max(vF, vH);

            // Load previous column for next segment's diagonal
            vH = pvHLoad[j];
        }

        // Lazy F correction
        for (int k = 0; k < VECLEN; k++) {
            vF.push(0);
            for (int j = 0; j < segLen; j++) {
                vH = pvHStore[j];
                vH = max(vH, vF);
                vMaxColumn = max(vMaxColumn, vH);
                pvHStore[j] = vH;
                vH = saturating_subu(vH, vGapO);
                vF = saturating_subu(vF, vGapE);
                // Early termination
                if (equal(max(vector<int16, VECLEN>(vF - vH), vZero), vZero)) goto end_lazy;
            }
        }
end_lazy:
        vMaxScore = max(vMaxScore, vMaxColumn);
    }

    return reduce_max(vMaxScore);
}

void sw_striped_kernel(
    input_stream_int32 *in_seq1,
    input_stream_int32 *in_seq2,
    output_stream_int32 *out_score
) {
    static int8_t query[MAX_SEQ_LEN];
    static int8_t ref[MAX_SEQ_LEN];

    int32_t pair_count = readincr(in_seq1);

    for (int pair = 0; pair < pair_count; pair++) {
        int32_t lens = readincr(in_seq1);
        int32_t qlen = (lens >> 16) & 0xFFFF;
        int32_t rlen = lens & 0xFFFF;

        // Read query — 2-bit encoding: 16 bases per int32
        for (int i = 0; i < (qlen + PACKETLEN - 1) / PACKETLEN; i++) {
            int32_t packed = readincr(in_seq1);
            int remaining = ((i * PACKETLEN + PACKETLEN) <= qlen) ? PACKETLEN : (qlen - i * PACKETLEN);
            for (int j = 0; j < remaining; j++)
                query[i * PACKETLEN + j] = (packed >> (j * 2)) & 0x03;
        }

        // Read ref — 2-bit encoding: 16 bases per int32
        int32_t sync = readincr(in_seq2);  // sync word (discarded)
        for (int i = 0; i < (rlen + PACKETLEN - 1) / PACKETLEN; i++) {
            int32_t packed = readincr(in_seq2);
            int remaining = ((i * PACKETLEN + PACKETLEN) <= rlen) ? PACKETLEN : (rlen - i * PACKETLEN);
            for (int j = 0; j < remaining; j++)
                ref[i * PACKETLEN + j] = (packed >> (j * 2)) & 0x03;
        }

        // Build profile and align
        vector<int16, VECLEN> vProfile[AA_LEN * SEG_LEN];
        build_profile(query, qlen, BLOSUM62_MAT, vProfile);
        int32_t score = sw_striped_align(ref, rlen, qlen, vProfile);

        writeincr(out_score, score);
    }
}
