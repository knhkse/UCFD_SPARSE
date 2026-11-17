#pragma once

#include "ucfdpcimpl.h"
#include "ucfdmatimpl.h"

typedef struct {
    UCFDInt bn;
    UCFDInt block;
    UCFDReal *diagvalues;
} Precon_BLUSGS;

typedef struct {
    Precon_BLUSGS base;
    UCFDInt ncolors;
    UCFDInt *icolors;
} Precon_PBLUSGS;

#if defined(__CUDACC__)
typedef struct {
    Precon_BLUSGS   base;
    cudaStream_t    stream;
    UCFDInt         ncolors;
    UCFDInt         *icolors;
    UCFDInt         *d_rowptr;
    UCFDInt         *d_colidx;
    UCFDInt         *d_diagslots;
} Precon_CUBLUSGS;
#endif

#if defined(__cplusplus)
extern "C" {
#endif

UCFD_INTERN ucfd_status_t BLUSGSPreconPrepare(Precon);
UCFD_INTERN ucfd_status_t BLUSGSPreconDestroy(Precon);

#if defined(__cplusplus)
}
#endif