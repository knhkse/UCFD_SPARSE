#pragma once

#include "ucfdsolverimpl.h"

typedef struct {
    UCFDInt restart;    // Restart number
    UCFDReal *H, *V, *y, *w, *sn, *cs, *htmp, *r;      // Working arrays
} Solver_GMRES;


#if defined(USE_CUDA)
typedef struct {
    cublasHandle_t handle;
    UCFDInt restart;
    UCFDInt ldv;                                // Padded leading dimension of V (>=n, even)
    UCFDReal *d_V, *d_y, *d_w, *d_r, *d_proj;   // Device arrays
    UCFDReal *H, *g, *y, *proj_host;            // Host arrays
} Solver_CUDAGMRES;
#endif