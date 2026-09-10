#pragma once

#include "ucfdtypes.h"


#if defined(UCFD_FLOAT32)
    #ifndef BETAST
        #define BETAST 0.09f
    #endif
    #ifndef GAMMA
        #define GAMMA 1.4f
    #endif
    #ifndef PMIN
        #define PMIN 1e-13f
    #endif
#else
    #ifndef BETAST
        #define BETAST 0.09
    #endif
    #ifndef GAMMA
        #define GAMMA 1.4
    #endif
    #ifndef PMIN
        #define PMIN 1e-13
    #endif
#endif

#if defined(__cplusplus)
extern "C" {
#endif

typedef enum {
    RANS_KWSST   = 22,
    RANS_SA      = 23
} rans_t;

typedef void (*fluxfunc)(UCFDInt, UCFDInt, UCFDInt, UCFDReal*, UCFDReal*, UCFDReal*);
typedef void (*srcjacobian)(UCFDInt, UCFDInt, UCFDReal*, UCFDReal*, UCFDReal*);

void ns_flux_container(UCFDInt nfvars, UCFDInt nturbvars, UCFDInt ndims, UCFDReal *u, UCFDReal *nf, UCFDReal *f);
void rans_flux_container(UCFDInt nfvars, UCFDInt nturbvars, UCFDInt ndims, UCFDReal *u, UCFDReal *nf, UCFDReal *f);

void kwsst_src_jacobian(UCFDInt nvars, UCFDInt nturbvars, UCFDReal *uf, UCFDReal *A, UCFDReal *dsrc);
void sa_src_jacobian(UCFDInt nvars, UCFDInt nturbvars, UCFDReal *uf, UCFDReal *A, UCFDReal *dsrc);

#if defined(__CUDACC__)
__device__ void cuda_kwsst_src_jacobian(UCFDInt nvars, UCFDInt nturbvars, UCFDReal *uf, UCFDReal *A, UCFDReal *dsrc);
__device__ void cuda_sa_src_jacobian(UCFDInt nvars, UCFDInt nturbvars, UCFDReal *uf, UCFDReal *A, UCFDReal *dsrc);
#endif


#if defined(__cplusplus)
}
#endif
