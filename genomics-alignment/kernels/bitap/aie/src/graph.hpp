#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <adf.h>
#include "kernels.hpp"

using namespace adf;

class bitapGraph : public graph {
    private:
        kernel ba;

    public:
        input_plio p_seq1;
        input_plio p_seq2;
        output_plio p_score;

        bitapGraph() {
            ba = kernel::create(bitap_kernel);
            source(ba) = "bitap_kernel.cc";

            p_seq1  = input_plio::create("StreamIn0", plio_32_bits, "data/seq1.txt");
            p_seq2  = input_plio::create("StreamIn1", plio_32_bits, "data/seq2.txt");
            p_score = output_plio::create("StreamOut0", plio_32_bits, "output_score.txt");

            connect<stream>(p_seq1.out[0], ba.in[0]);
            connect<stream>(p_seq2.out[0], ba.in[1]);
            connect<stream>(ba.out[0], p_score.in[0]);

            runtime<ratio>(ba) = 0.9;
        };
};

#endif
