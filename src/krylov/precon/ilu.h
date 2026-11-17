#pragma once

#include "ucfdpcimpl.h"
#include "ucfdmatimpl.h"

typedef struct {
    UCFDInt     *iw;
    UCFDInt     n;
} Precon_ILU;

typedef struct {
    UCFDInt     *iw;
    UCFDInt     bn;
    UCFDInt     block;
} Precon_BILU;


typedef struct {
    Precon_BILU base;
    UCFDInt     ncolors;
    UCFDInt     *icolors;
} Precon_PBILU;

#if defined(__CUDACC__)
typedef struct {
    Precon_BILU     base;
    cudaStream_t    stream;
    UCFDInt         ncolors;
    UCFDInt         *icolors;
    UCFDInt         *d_rowptr;
    UCFDInt         *d_colidx;
    UCFDInt         *d_diagslots;
} Precon_CUBILU;
#endif

#if defined(__cplusplus)
extern "C" {
#endif

/* Destroy inner structure */
UCFD_INTERN ucfd_status_t ILUPreconDestroy(Precon precon);

#if defined(__cplusplus)
}
#endif
