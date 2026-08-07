#ifndef COMPUTE_HH
#define COMPUTE_HH

// Compute kernel router — selected at build time via KERNEL_CONFIG define
// All compute kernels implement:
//   void compute_init(int32_t* batch, int batch_size, SwProfile& prof);
//   int32_t compute_score(int32_t* batch, int batch_size, int item_index, int row, int col,
//                         const SwProfile& prof);

#if defined(USE_WFA)
#include "compute_wfa.hh"
#elif defined(USE_BITAP)
#include "compute_bitap.hh"
#elif defined(USE_STRIPED_SWG)
#include "compute_striped_swg.hh"
#elif defined(USE_INTERSEQ_SWG)
#include "compute_swg_interseq.hh"
#elif defined(USE_SCALAR_SWG)
#include "compute_swg_scalar.hh"
#else
#include "compute_dummy.hh"
#endif

#endif
