#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#include "sparseformat.h"


__global__ static void
_SpMV_CUDACSR(UCFDReal alpha, SpMat A, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    UCFDInt i, rs, ed, j;
    UCFDReal s;

    i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < A->n)
    {
        rs = A->rowptr[i];
        ed = A->rowptr[i+1];
        s = 0.0;
        for (j=rs; j<ed; j++) {
            s += A->values[j] * x[A->colidx[j]];
        }
        y[i] = alpha*s + beta*y[i];
    }
}


static inline ucfd_status_t
SpMV_CUDACSR(UCFDReal alpha, SpMat A, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    UCFDInt bpg = (A->n + TPB - 1)/TPB;
    _SpMV_CUDACSR<<<bpg, TPB>>>(alpha, A, x, beta, y);
    UCFDFunctionReturn(UCFD_SUCCESS);
}



ucfd_status_t MatCreateCUDACSR(SpMat *mat, UCFDInt n, UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values)
{
    UCFDCall(UCFDMatInit(mat));

}


















