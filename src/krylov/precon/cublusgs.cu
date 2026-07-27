#include <cuda_runtime.h>

#include "blusgs.h"
#include "inverse_cuda.cuh"


__global__ static void
CUDABLUSGSPreconPreparePerColor(UCFDInt bn, UCFDInt block,
                                UCFDInt __restrict__ *diagslots,
                                UCFDReal __restrict__ *diagvalues,
                                UCFDReal __restrict__ *values)
{
    UCFDInt idx, didx;
    UCFDInt blkdim = block*block;
    UCFDReal *diagblock;

    idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < bn)
    {
        didx = diagslots[idx];
        diagblock = &diagvalues[idx*blkdim];
        CUDACall(cudaMemcpy(diagblock, &values[didx*blkdim], blkdim*sizeof(UCFDReal),
                            cudaMemcpyDeviceToDevice));
        ludcmp(block, diagblock);
    }
}

static ucfd_status_t
CUDABLUSGSPreconPrepare(Precon precon)
{
    Precon_PBLUSGS *pblusgs = (Precon_PBLUSGS *)precon->data;
    UCFDInt bpg = (pblusgs->base.bn + TPB - 1)/TPB;

    CUDABLUSGSPreconPreparePerColor<<<bpg, TPB>>>(
        pblusgs->base.bn, pblusgs->base.block, precon->diagslots,
        pblusgs->base.diagvalues, precon->values
    );
    UCFDFunctionReturn(UCFD_SUCCESS);
}


/* Forward sweep : (D+L)x' = b -> x' = inv(D) * (b-Lx') */
template<UCFDInt block>
__global__ static void
CUDABLUSGSPreconLowerApply(UCFDInt nstart, UCFDInt nend,
                           UCFDInt __restrict__ *rowptr, UCFDInt __restrict__ *colidx,
                           UCFDReal __restrict__ *values, UCFDInt __restrict__ *diagslots,
                           UCFDReal __restrict__ *diagvalues, UCFDReal __restrict__ *b)
{
    UCFDInt idx, jdx, kdx, row, col;
    UCFDInt dd, st, cind;
    UCFDInt blkdim = block*block;
    UCFDInt range = nend - nstart;
    UCFDReal v, arr[block];

    idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < range)
    {
        idx += nstart;
        dd = diagslots[idx];
        st = rowptr[idx];

        // Initialize arr
        for (kdx=0; kdx<block; kdx++)
            arr[kdx] = b[kdx + idx*block];

        // arr := b - Lx'
        for (jdx=st; jdx<dd; jdx++)
        {
            cind = colidx[jdx];
            for (row = 0; row < block; row++)
            {
                v = 0.0;
                for (col = 0; col < block; col++)
                    v += values[col + row * block + jdx * blkdim] * b[col + cind * block];
                arr[row] -= v;
            }
        }

        // x' := inv(D) * (b - Lx') = inv(D)*arr
        lusub(block, &diagvalues[idx*blkdim], arr);
        for (kdx=0; kdx<block; kdx++)
            b[kdx + idx*block] = arr[kdx];
    }
}


/* Backward sweep : (D+U)x = Dx' -> x = x' - inv(D) * Ux */
template<UCFDInt block>
__global__ static void
CUDABLUSGSPreconUpperApply(UCFDInt nstart, UCFDInt nend,
                           UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values,
                           UCFDInt *diagslots, UCFDReal *diagvalues, UCFDReal *b)
{
    UCFDInt idx, jdx, kdx, row, col;
    UCFDInt dd, ed, cind;
    UCFDInt blkdim = block*block;
    UCFDInt range = nend - nstart;
    UCFDReal v, arr[block];

    idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < range)
    {
        idx += nstart;
        dd = diagslots[idx];
        ed = rowptr[idx+1];

        // Initialize
        for (kdx = 0; kdx < block; kdx++)
            arr[kdx] = 0.0;
        
        // arr := Ux
        for (jdx = dd + 1; jdx < ed; jdx++)
        {
            cind = colidx[jdx];
            for (row = 0; row < block; row++)
            {
                v = 0.0;
                for (col = 0; col < block; col++)
                    v += values[col + row * block + jdx * blkdim] * b[col + cind * block];
                arr[row] += v;
            }
        }

        // arr := inv(D) * Ux
        lusub(block, &diagvalues[idx*blkdim], arr);

        // b := inv(D) * Ux
        for (kdx = 0; kdx < block; kdx++)
            b[kdx + idx * block] -= arr[kdx];
    }
}

static ucfd_status_t
CUDABLUSGSPreconApply(Precon precon, UCFDReal *b)
{
    Precon_PBLUSGS *pblusgs = (Precon_PBLUSGS *)precon->data;
    UCFDInt i;
    UCFDInt bpg = (pblusgs->base.bn + TPB - 1)/TPB;

    switch (pblusgs->base.block) {
        case 1:
            // UCFDWarning("Single block size(block=1)::Use LU-SGS preconditioner")
            for (i=0; i<pblusgs->ncolors; i++)
                CUDABLUSGSPreconLowerApply<1><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            for (i=pblusgs->ncolors-1; i>=0; i--)
                CUDABLUSGSPreconUpperApply<1><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            break;
        case 2:
            for (i=0; i<pblusgs->ncolors; i++)
                CUDABLUSGSPreconLowerApply<2><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            for (i=pblusgs->ncolors-1; i>=0; i--)
                CUDABLUSGSPreconUpperApply<2><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            break;
        case 3:
            for (i=0; i<pblusgs->ncolors; i++)
                CUDABLUSGSPreconLowerApply<3><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            for (i=pblusgs->ncolors-1; i>=0; i--)
                CUDABLUSGSPreconUpperApply<3><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            break;
        case 4:
            for (i=0; i<pblusgs->ncolors; i++)
                CUDABLUSGSPreconLowerApply<4><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            for (i=pblusgs->ncolors-1; i>=0; i--)
                CUDABLUSGSPreconUpperApply<4><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            break;
        case 5:
            for (i=0; i<pblusgs->ncolors; i++)
                CUDABLUSGSPreconLowerApply<5><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            for (i=pblusgs->ncolors-1; i>=0; i--)
                CUDABLUSGSPreconUpperApply<5><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            break;
        case 6:
            for (i=0; i<pblusgs->ncolors; i++)
                CUDABLUSGSPreconLowerApply<6><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            for (i=pblusgs->ncolors-1; i>=0; i--)
                CUDABLUSGSPreconUpperApply<6><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            break;
        case 7:
            for (i=0; i<pblusgs->ncolors; i++)
                CUDABLUSGSPreconLowerApply<7><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            for (i=pblusgs->ncolors-1; i>=0; i--)
                CUDABLUSGSPreconUpperApply<7><<<bpg, TPB>>>(
                    pblusgs->icolors[i], pblusgs->icolors[i+1],
                    precon->rowptr, precon->colidx, precon->values,
                    precon->diagslots, pblusgs->base.diagvalues, b
                );
            break;
        default: fprintf(stderr, "Unsupported block size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t CUDABLUSGSPreconDestroy(Precon precon)
{
    if (!precon) UCFDFunctionReturn(UCFD_SUCCESS);
    Precon_BLUSGS *blusgs = (Precon_BLUSGS *)precon->data;
    cudaFree(blusgs->diagvalues);
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDPreconSetCUDABLUSGS(Precon *precon, UCFDInt bn, UCFDInt block, UCFDInt ncolors, UCFDInt *icolors)
{
    Precon pc = *precon;

    CheckCUDAPointer(pc->rowptr);
    CheckCUDAPointer(pc->colidx);
    CheckCUDAPointer(pc->values);
    CheckCUDAPointer(pc->diagslots);

    Precon_PBLUSGS *cudablusgs = (Precon_PBLUSGS *)calloc(1, sizeof(*cudablusgs));
    UCFDCheckNull(cudablusgs, "CUDA BLUSGS precon allocation failed\n");

    cudablusgs->base.bn         = bn;
    cudablusgs->base.block      = block;
    cudablusgs->ncolors         = ncolors;
    cudablusgs->icolors         = icolors;
    CUDACall(cudaMalloc((void**)&cudablusgs->base.diagvalues, (size_t)(bn*block*block)*sizeof(UCFDReal)));

    pc->type_name       = CUDABLUSGS;
    pc->data            = cudablusgs;
    pc->ops->prepare    = CUDABLUSGSPreconPrepare;
    pc->ops->apply      = CUDABLUSGSPreconApply;
    pc->ops->destroy    = CUDABLUSGSPreconDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
