#include <cuda_runtime.h>

#include "blusgs.h"
#include "inverse_cuda.cuh"


__global__ static void
CUDABLUSGSPreconPreparePerColor(UCFDInt bn, UCFDInt block,
                                const UCFDInt *__restrict__ diagslots,
                                const UCFDReal *__restrict__ values,
                                UCFDReal *__restrict__ diagvalues)
{
    const UCFDInt idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= bn) return;
    
    const UCFDInt blkdim = block*block;
    const UCFDInt didx = diagslots[idx];
    UCFDReal *diagblock = &diagvalues[idx*blkdim];

    for (UCFDInt i=0; i<blkdim; ++i)
        diagblock[i] = values[didx*blkdim + i];

    ludcmp(block, diagblock);
}

static ucfd_status_t
CUDABLUSGSPreconPrepare(Precon precon)
{
    Precon_CUBLUSGS *cublusgs = (Precon_CUBLUSGS *)precon->data;
    const UCFDInt bpg = (cublusgs->base.bn + TPB - 1)/TPB;

    CUDABLUSGSPreconPreparePerColor<<<bpg, TPB, 0, cublusgs->stream>>>(
        cublusgs->base.bn, cublusgs->base.block, cublusgs->d_diagslots,
        precon->values, cublusgs->base.diagvalues
    );
    UCFDFunctionReturn(UCFD_SUCCESS);
}


/* Forward sweep : (D+L)x' = b -> x' = inv(D) * (b-Lx') */
template<UCFDInt block>
__global__ static void
CUDABLUSGSPreconLowerApply(const UCFDInt interval, const UCFDInt nstart,
                           const UCFDInt *__restrict__ rowptr, 
                           const UCFDInt *__restrict__ colidx,
                           const UCFDReal *__restrict__ values,
                           const UCFDInt *__restrict__ diagslots,
                           UCFDReal *__restrict__ diagvalues,
                           UCFDReal *__restrict__ b)
{
    const UCFDInt _idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (_idx >= interval) return;
    
    const UCFDInt blkdim = block*block;
    UCFDInt jdx, kdx, row, col, cind;
    UCFDReal v, arr[block];

    const UCFDInt idx = _idx + nstart;
    const UCFDInt dd = diagslots[idx];
    const UCFDInt st = rowptr[idx];

    // Initialize arr
    for (kdx=0; kdx<block; ++kdx)
        arr[kdx] = b[kdx + idx*block];

    // arr := b - Lx'
    for (jdx=st; jdx<dd; ++jdx)
    {
        cind = colidx[jdx];
        for (row = 0; row < block; ++row)
        {
            v = 0.0;
            for (col = 0; col < block; ++col)
                v += values[col + row * block + jdx * blkdim] * b[col + cind * block];
            arr[row] -= v;
        }
    }

    // x' := inv(D) * (b - Lx') = inv(D)*arr
    lusub(block, &diagvalues[idx*blkdim], arr);
    for (kdx=0; kdx<block; ++kdx)
        b[kdx + idx*block] = arr[kdx];
}


/* Backward sweep : (D+U)x = Dx' -> x = x' - inv(D) * Ux */
template<UCFDInt block>
__global__ static void
CUDABLUSGSPreconUpperApply(UCFDInt interval, UCFDInt nstart,
                           const UCFDInt *__restrict__ rowptr, 
                           const UCFDInt *__restrict__ colidx,
                           const UCFDReal *__restrict__ values,
                           const UCFDInt *__restrict__ diagslots,
                           UCFDReal *__restrict__ diagvalues,
                           UCFDReal *__restrict__ b)
{
    const UCFDInt _idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (_idx >= interval) return;
    
    const UCFDInt blkdim = block*block;
    UCFDInt jdx, kdx, row, col, cind;
    UCFDReal v, arr[block];

    const UCFDInt idx = _idx + nstart;
    const UCFDInt dd = diagslots[idx];
    const UCFDInt ed = rowptr[idx+1];

    // Initialize
    for (kdx = 0; kdx < block; ++kdx)
        arr[kdx] = 0.0;
    
    // arr := Ux
    for (jdx = dd + 1; jdx < ed; ++jdx)
    {
        cind = colidx[jdx];
        for (row = 0; row < block; ++row)
        {
            v = 0.0;
            for (col = 0; col < block; ++col)
                v += values[col + row * block + jdx * blkdim] * b[col + cind * block];
            arr[row] += v;
        }
    }

    // arr := inv(D) * Ux
    lusub(block, &diagvalues[idx*blkdim], arr);

    // b := inv(D) * Ux
    for (kdx = 0; kdx < block; ++kdx)
        b[kdx + idx * block] -= arr[kdx];
}

static ucfd_status_t
CUDABLUSGSPreconApply(Precon precon, UCFDReal *b)
{
    Precon_CUBLUSGS *cublusgs = (Precon_CUBLUSGS *)precon->data;
    UCFDInt i, nstart, interval, bpg;
    const UCFDInt ncolors = cublusgs->ncolors;

    switch (cublusgs->base.block) {
        case 1:
            for (i=0; i<ncolors; ++i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconLowerApply<1><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconUpperApply<1><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            break;
        case 2:
            for (i=0; i<ncolors; ++i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconLowerApply<2><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconUpperApply<2><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            break;
        case 3:
            for (i=0; i<ncolors; ++i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconLowerApply<3><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconUpperApply<3><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            break;
        case 4:
            for (i=0; i<ncolors; ++i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconLowerApply<4><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconUpperApply<4><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            break;
        case 5:
            for (i=0; i<ncolors; ++i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconLowerApply<5><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconUpperApply<5><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            break;
        case 6:
            for (i=0; i<ncolors; ++i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconLowerApply<6><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconUpperApply<6><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            break;
        case 7:
            for (i=0; i<ncolors; ++i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconLowerApply<7><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = cublusgs->icolors[i];
                interval = cublusgs->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABLUSGSPreconUpperApply<7><<<bpg, TPB, 0, cublusgs->stream>>>(
                    interval, nstart,
                    cublusgs->d_rowptr, cublusgs->d_colidx, precon->values,
                    cublusgs->d_diagslots, cublusgs->base.diagvalues, b
                );
            }
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
UCFDPreconSetCUDABLUSGS(Precon *precon, SpMat mat, UCFDInt bn, UCFDInt block, UCFDInt ncolors, UCFDInt *icolors)
{
    UCFDCheckNull(*precon, "Preconditioner must be initialized\n");
    Precon pc = *precon;
    Precon_CUBLUSGS *cublusgs = (Precon_CUBLUSGS *)calloc(1, sizeof(*cublusgs));
    UCFDCheckNull(cublusgs, "CUDA BLUSGS precon allocation failed\n");

    cublusgs->base.bn         = bn;
    cublusgs->base.block      = block;
    cublusgs->ncolors         = ncolors;
    cublusgs->icolors         = icolors;
    pc->type_name             = CUDABLUSGS;

    UCFDCall(UCFDCopyStream(mat, &cublusgs->stream));
    CUDACall(cudaMalloc(&cublusgs->d_rowptr, (bn+1)*sizeof(UCFDInt)));
    CUDACall(cudaMalloc(&cublusgs->d_colidx, pc->nnz*sizeof(UCFDInt)));
    CUDACall(cudaMalloc(&cublusgs->d_diagslots, bn*sizeof(UCFDInt)));
    CUDACall(cudaMalloc(&pc->values, pc->nnz*block*block*sizeof(UCFDReal)));
    CUDACall(cudaMalloc(&cublusgs->base.diagvalues, (size_t)(bn*block*block)*sizeof(UCFDReal)));
    CUDACall(cudaMemcpy(cublusgs->d_rowptr, pc->rowptr, (bn+1)*sizeof(UCFDInt), cudaMemcpyHostToDevice));
    CUDACall(cudaMemcpy(cublusgs->d_colidx, pc->colidx, pc->nnz*sizeof(UCFDInt), cudaMemcpyHostToDevice));
    CUDACall(cudaMemcpy(cublusgs->d_diagslots, pc->diagslots, bn*sizeof(UCFDInt), cudaMemcpyHostToDevice));

    pc->data            = cublusgs;
    pc->ops->prepare    = CUDABLUSGSPreconPrepare;
    pc->ops->apply      = CUDABLUSGSPreconApply;
    pc->ops->destroy    = CUDABLUSGSPreconDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
