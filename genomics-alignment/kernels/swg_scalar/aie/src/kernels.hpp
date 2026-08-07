#ifndef __KERNELS_H__
#define __KERNELS_H__

#include <adf/stream/types.h>

void swg_kernel(input_stream_int32 *in_seq1, input_stream_int32 *in_seq2, output_stream_int32 *out_score);

#endif
