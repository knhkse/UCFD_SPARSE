#pragma once

#include "ucfdsolverimpl.h"

typedef struct {
    UCFDInt     n;
    UCFDInt     restart;    // Restart number
    UCFDReal    *H, *V, *y, *w, *sn, *cs, *htmp, *r;      // Working arrays
} Solver_GMRES;

#if defined(__CUDACC__)
typedef struct {
    UCFDInt         n;
    cublasHandle_t  handle;
    cudaStream_t    stream;
    UCFDInt         restart;
    UCFDInt         ldv;                    /* Padded leading dimension of V (>=n, even) */
    UCFDReal        *d_V, *d_y, *d_proj;    /* Device arrays */
    UCFDReal        *H, *g;
    UCFDReal        *y, *proj_host;         /* Pinned */
} Solver_CUDAGMRES;
#endif
