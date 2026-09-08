#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/fill.h>

#include "ilu.h"
#include "inverse_cuda.cuh"


template<UCFDInt block>
__global__ static void
CUDABILUPreconPreparePerColor(const UCFDInt interval, const UCFDInt nstart,
                              const UCFDInt *__restrict__ rowptr,
                              const UCFDInt *__restrict__ colidx,
                              const UCFDInt *__restrict__ diagslots,
                              UCFDReal *__restrict__ values,
                              UCFDInt *__restrict__ iw)
{
    const UCFDInt _idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (_idx >= interval) return;

    const UCFDInt blkdim = block*block;
    UCFDInt kdx, row, col, ele;
    UCFDInt ck, kk, kst, ked, jj, iwj;
    UCFDReal v, Aik[block][block];

    const UCFDInt idx = _idx + nstart;
    const UCFDInt st = rowptr[idx];
    const UCFDInt jed = rowptr[idx+1];
    const UCFDInt ed = diagslots[idx];

    for (kdx=st; kdx<ed; ++kdx)
    {
        ck = colidx[kdx];
        kk = diagslots[ck];
        kst = rowptr[ck];
        ked = rowptr[ck+1];

        // A[i,k] := A[i,k] @ inv(A[k,k])
        lusubmattrans(block, &(values[kk*blkdim]), &(values[kdx*blkdim]));
        // memcpy(Aik, &values[kdx*blkdim], sizeof(double)*block);
        for (row=0; row<block; ++row) {
            for (col=0; col<block; ++col)
                Aik[row][col] = values[kdx*blkdim+row*block+col];
        }

        // Prepare iw
        for (jj=kst; jj<ked; ++jj) iw[colidx[jj]] = jj;

        // j iteration
        for (jj=kdx+1; jj<jed; ++jj) {
            iwj = iw[colidx[jj]];

            if (iwj != -1) {
                // values[jj] -= Aik * values[iwj]
                for (row=0; row<block; ++row) {
                    for (col=0; col<block; ++col) {
                        v = 0.0;
                        for (ele=0; ele<block; ++ele)
                            v += Aik[row][ele] * values[iwj*blkdim+ele*block+col];
                        values[jj*blkdim+row*block+col] -= v;
                    }
                }
            }
        }
        // Clean iw
        for (jj=kst; jj<ked; ++jj) iw[colidx[jj]] = -1;
    }
    // LU decomposition of current row diagonal matrix
    ludcmp(block, &(values[ed*blkdim]));
}

static ucfd_status_t
CUDABILUPreconPrepare(Precon precon)
{
    Precon_PBILU *pbilu = (Precon_PBILU *)precon->data;
    UCFDInt i, nstart, interval, bpg;
    const UCFDInt ncolors = pbilu->ncolors;

    switch (pbilu->base.block) {
        case 1:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconPreparePerColor<1><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->diagslots, precon->values, pbilu->base.iw 
                );
            }
            break;
        case 2:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconPreparePerColor<2><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->diagslots, precon->values, pbilu->base.iw 
                );
            }
            break;
        case 3:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconPreparePerColor<3><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->diagslots, precon->values, pbilu->base.iw 
                );
            }
            break;
        case 4:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconPreparePerColor<4><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->diagslots, precon->values, pbilu->base.iw 
                );
            }
            break;
        case 5:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconPreparePerColor<5><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->diagslots, precon->values, pbilu->base.iw 
                );
            }
            break;
        case 6:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconPreparePerColor<6><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->diagslots, precon->values, pbilu->base.iw 
                );
            }
            break;
        case 7:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconPreparePerColor<7><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->diagslots, precon->values, pbilu->base.iw 
                );
            }
            break;
        default: fprintf(stderr, "Unsupported block size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

template<UCFDInt block>
__global__ static void
CUDABILUPreconLowerApply(const UCFDInt interval, const UCFDInt nstart,
                         const UCFDInt *__restrict__ rowptr,
                         const UCFDInt *__restrict__ colidx,
                         const UCFDReal *__restrict__ values,
                         const UCFDInt *__restrict__ diagslots,
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
    for (kdx = 0; kdx < block; ++kdx)
        arr[kdx] = b[kdx + idx * block];

    for (jdx = st; jdx < dd; ++jdx)
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
    for (kdx = 0; kdx < block; ++kdx)
        b[kdx + idx * block] = arr[kdx];
}

template<UCFDInt block>
__global__ static void
CUDABILUPreconUpperApply(const UCFDInt interval, const UCFDInt nstart,
                         const UCFDInt *__restrict__ rowptr,
                         const UCFDInt *__restrict__ colidx,
                         const UCFDReal *__restrict__ values,
                         const UCFDInt *__restrict__ diagslots,
                         UCFDReal *__restrict__ b)
{
    const UCFDInt _idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (_idx >= interval) return;

    const UCFDInt blkdim = block*block;
    UCFDInt jdx, kdx, row, col, cind;
    UCFDReal v, arr[block];

    const UCFDInt idx = _idx + nstart;
    const UCFDInt dd = diagslots[idx];
    const UCFDInt ed = rowptr[idx + 1];

    // Initialize
    for (kdx = 0; kdx < block; ++kdx)
        arr[kdx] = b[kdx + idx * block];

    for (jdx = dd + 1; jdx < ed; ++jdx)
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
    // LU substitution for vector
    lusub(block, &(values[dd*blkdim]), arr);
    for (row=0; row<block; ++row)
        b[idx*block+row] = arr[row];
}

static ucfd_status_t
CUDABILUPreconApply(Precon precon, UCFDReal *b)
{
    Precon_PBILU *pbilu = (Precon_PBILU *)precon->data;
    UCFDInt i, nstart, interval, bpg;
    const UCFDInt ncolors = pbilu->ncolors;


    switch (pbilu->base.block) {
        case 1:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconLowerApply<1><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconUpperApply<1><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            break;
        case 2:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconLowerApply<2><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconUpperApply<2><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            break;
        case 3:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconLowerApply<3><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconUpperApply<3><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            break;
        case 4:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconLowerApply<4><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconUpperApply<4><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            break;
        case 5:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconLowerApply<5><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconUpperApply<5><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            break;
        case 6:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconLowerApply<6><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconUpperApply<6><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            break;
        case 7:
            for (i=0; i<ncolors; ++i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconLowerApply<7><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            for (i=ncolors-1; i>=0; --i) {
                nstart = pbilu->icolors[i];
                interval = pbilu->icolors[i+1] - nstart;
                bpg = (interval + TPB - 1)/TPB;
                CUDABILUPreconUpperApply<7><<<bpg, TPB>>>(
                    interval, nstart, precon->rowptr,
                    precon->colidx, precon->values, precon->diagslots, b
                );
            }
            break;
        default:
            fprintf(stderr, "Unsupported block size\n");
            UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t CUDABILUPreconDestroy(Precon precon)
{
    if (!precon) UCFDFunctionReturn(UCFD_SUCCESS);
    Precon_BILU *bilu = (Precon_BILU *)precon->data;
    cudaFree(bilu->iw);
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDPreconSetCUDABILU(Precon *precon, UCFDInt bn, UCFDInt block, UCFDInt ncolors, UCFDInt *icolors)
{
    Precon pc = *precon;

    CheckCUDAPointer(pc->rowptr);
    CheckCUDAPointer(pc->colidx);
    CheckCUDAPointer(pc->values);
    CheckCUDAPointer(pc->diagslots);

    Precon_PBILU *cudabilu = (Precon_PBILU *)calloc(1, sizeof(*cudabilu));
    UCFDCheckNull(cudabilu, "CUDA BILU precon allocation failed\n");

    /* Initialize working array */
    CUDACall(cudaMalloc((void**)&((Precon_BILU *)cudabilu)->iw, (size_t)bn*sizeof(UCFDInt)));
    thrust::device_ptr<UCFDInt> dev_ptr(((Precon_BILU *)cudabilu)->iw);
    thrust::fill(dev_ptr, dev_ptr+bn, -1);

    ((Precon_BILU *)cudabilu)->bn          = bn;
    ((Precon_BILU *)cudabilu)->block       = block;
    cudabilu->ncolors                      = ncolors;
    cudabilu->icolors                      = icolors;

    pc->type_name                       = CUDABILU;
    CUDACall(cudaMalloc(&pc->values, pc->nnz*block*block*sizeof(UCFDReal)));
    pc->data                            = cudabilu;
    pc->ops->prepare                    = CUDABILUPreconPrepare;
    pc->ops->apply                      = CUDABILUPreconApply;
    pc->ops->destroy                    = CUDABILUPreconDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
