#ifndef __KERNELS_HPP__
#define __KERNELS_HPP__

#include <adf.h>
#include "aie_api/aie.hpp"
#include <aie_api/aie_adf.hpp>

void bitap_kernel(
    input_stream_int32 *in_seq1,
    input_stream_int32 *in_seq2,
    output_stream_int32 *out_score
);

#endif
