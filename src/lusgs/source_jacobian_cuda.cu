#include "flux.h"

__device__ void
cuda_kwsst_src_jacobian(UCFDInt nvars, UCFDInt nturbvars, UCFDReal *uf, UCFDReal *A, UCFDReal *dsrc)
{
    const UCFDReal k = uf[nvars-2]/uf[0];
    A[0] += dsrc[nvars-2];
    A[1] += max(BETAST*k, 0.0);
    A[1 + nturbvars] = dsrc[nvars-1];
}

__device__ void
cuda_sa_src_jacobian(UCFDInt nvars, UCFDInt nturbvars, UCFDReal *uf, UCFDReal *A, UCFDReal *dsrc)
{
    A[0] += dsrc[nvars-1];
}
