#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#include "sparsemat.h"


__global__ static void
_SpMV_CUDACSR(UCFDReal alpha, UCFDInt n,
              UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values,
              UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    UCFDInt idx, jdx, rs, ed;
    UCFDReal s;

    idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < n)
    {
        rs = rowptr[idx];
        ed = rowptr[idx+1];
        s = 0.0;
        for (jdx=rs; jdx<ed; jdx++) {
            s += values[jdx] * x[colidx[jdx]];
        }
        y[idx] = alpha*s + beta*y[idx];
    }
}


static inline ucfd_status_t
SpMV_CUDACSR(UCFDReal alpha, SpMat mat, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
#if defined(DEBUG)
    CheckCUDAPointer(x);
    CheckCUDAPointer(y);
#endif
    BaseCSR *A = (BaseCSR *)mat->data;
    UCFDInt bpg = (A->n + TPB - 1)/TPB;
    _SpMV_CUDACSR<<<bpg, TPB>>>(alpha, A->n, A->rowptr, A->colidx, A->values, x, beta, y);
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDMatCreateCUDACSR(SpMat *mat, UCFDInt n, UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values)
{
    UCFDCall(UCFDMatInit(mat));
    SpMat m = *mat;
    m->type_name = CSRCUDA;

    BaseCSR *csr = (BaseCSR *)calloc(1, sizeof(*csr));
    UCFDCheckNull(csr, "CSR matrix creation failed\n");

    /* Check if input arrays are allocated in device side */
    CheckCUDAPointer(rowptr);
    CheckCUDAPointer(colidx);
    CheckCUDAPointer(values);

    csr->n              = n;
    csr->rowptr         = rowptr;
    csr->colidx         = colidx;
    csr->values         = values;
    m->data             = csr;
    m->ops->spmv        = SpMV_CUDACSR;
    m->ops->destroy     = UCFDEmptyKernel;

    UCFDFunctionReturn(UCFD_SUCCESS);
}


template<UCFDInt blk>
__global__ static void
_SpMV_CUDABSR(UCFDReal alpha, UCFDInt bn,
              UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values,
              UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    UCFDInt idx, jdx, kdx, row, col, st, ed;
    UCFDReal v, mat[blk*blk], arr[blk], xprt[blk];
    UCFDInt blk2 = blk*blk;

    idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx < bn)
    {
        /* Initialize */
        for (kdx=0; kdx<blk; kdx++) arr[kdx] = 0.0;
        st = rowptr[idx];
        ed = rowptr[idx+1];

        for (jdx=st; jdx<ed; jdx++)
        {
            kdx = colidx[jdx];

            for (row=0; row<blk; row++) {
                for (col=0; col<blk; col++) {
                    mat[row*blk+col] = values[jdx*blk2 + row*blk + col];
                }
                xprt[row] = x[kdx*blk+row];
            }
            // blockmv
            for (row=0; row<blk; row++) {
                v = 0.0;
                for (col=0; col<blk; col++)
                    v += mat[row*blk+col] * xprt[col];
                arr[row] += v;
            }
        }
        for (kdx=0; kdx<blk; kdx++)
            y[idx*blk+kdx] = alpha*arr[kdx] + beta*y[idx*blk+kdx];
    }
}


static ucfd_status_t
SpMV_CUDABSR(UCFDReal alpha, SpMat mat, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
#if defined(DEBUG)
    CheckCUDAPointer(x);
    CheckCUDAPointer(y);
#endif
    BaseCSR *A = (BaseCSR *)mat->data;
    BaseBSR *bsr = (BaseBSR *)mat->data;
    UCFDInt bpg = (bsr->bn + TPB - 1)/TPB;

    switch (bsr->block) {
        case 1: _SpMV_CUDABSR<1><<<bpg, TPB>>>(alpha, bsr->bn, A->rowptr, A->colidx, A->values, x, beta, y); break;
        case 2: _SpMV_CUDABSR<2><<<bpg, TPB>>>(alpha, bsr->bn, A->rowptr, A->colidx, A->values, x, beta, y); break;
        case 3: _SpMV_CUDABSR<3><<<bpg, TPB>>>(alpha, bsr->bn, A->rowptr, A->colidx, A->values, x, beta, y); break;
        case 4: _SpMV_CUDABSR<4><<<bpg, TPB>>>(alpha, bsr->bn, A->rowptr, A->colidx, A->values, x, beta, y); break;
        case 5: _SpMV_CUDABSR<5><<<bpg, TPB>>>(alpha, bsr->bn, A->rowptr, A->colidx, A->values, x, beta, y); break;
        case 6: _SpMV_CUDABSR<6><<<bpg, TPB>>>(alpha, bsr->bn, A->rowptr, A->colidx, A->values, x, beta, y); break;
        case 7: _SpMV_CUDABSR<7><<<bpg, TPB>>>(alpha, bsr->bn, A->rowptr, A->colidx, A->values, x, beta, y); break;
        default: fprintf(stderr, "Unsupported block size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDMatCreateCUDABSR(SpMat *mat, UCFDInt bn, UCFDInt blk, UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values)
{
    UCFDCall(UCFDMatInit(mat));
    SpMat m = *mat;
    m->type_name = BSRCUDA;

    BaseBSR *bsr = (BaseBSR *)calloc(1, sizeof(*bsr));
    UCFDCheckNull(bsr, "BSR matrix creation failed\n");

    bsr->bn             = bn;
    bsr->block          = blk;

    /* Check if input arrays are allocated in device side */
    CheckCUDAPointer(rowptr);
    CheckCUDAPointer(colidx);
    CheckCUDAPointer(values);

    ((BaseCSR *)bsr)->n         = (UCFDInt)(bn*blk);
    ((BaseCSR *)bsr)->rowptr    = rowptr;
    ((BaseCSR *)bsr)->colidx    = colidx;
    ((BaseCSR *)bsr)->values    = values;
    m->data                     = bsr;
    m->ops->spmv                = SpMV_CUDABSR;
    m->ops->destroy             = UCFDEmptyKernel;

    UCFDFunctionReturn(UCFD_SUCCESS);
}


#if defined(USE_CUSPARSE)
/**
 * Generic code for SpMV operation
 */
static ucfd_status_t SpMV_CUSPARSE(UCFDReal alpha, SpMat mat, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
#if defined(DEBUG)
    CheckCUDAPointer(x);
    CheckCUDAPointer(y);
#endif
    UCFDReal alpha_ = alpha, beta_ = beta;
    cuSPARSEContext *ctx = (cuSPARSEContext *)mat->data;
    CUSPARSECall(cusparseDnVecSetValues(ctx->vecX, x));
    CUSPARSECall(cusparseDnVecSetValues(ctx->vecY, y));
    CUSPARSECall(cusparseSpMV(
        ctx->handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha_,
        ctx->op, ctx->vecX, &beta_, ctx->vecY, CUSPARSE_REALTYPE,
        CUSPARSE_SPMV_ALG_DEFAULT, ctx->buffer
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t Destroy_CUSPARSE(SpMat mat)
{
    if (!mat) UCFDFunctionReturn(UCFD_SUCCESS);
    cuSPARSEContext *cumat = (cuSPARSEContext *)mat->data;

    CUSPARSECall(cusparseDestroy(cumat->handle));

    CUSPARSECall(cusparseDestroySpMat(cumat->op));
    CUSPARSECall(cusparseDestroyDnVec(cumat->vecX));
    CUSPARSECall(cusparseDestroyDnVec(cumat->vecY));
    CUDACall(cudaFree(cumat->buffer));
    CUDACall(cudaFree(cumat->tmp));

    UCFDFunctionReturn(UCFD_SUCCESS);
}

/**
 * Legacy code for BSR format SpMV operation
 * cusparseSpMV function supports BSR format after ver.13.0.1
 */
static ucfd_status_t
SpMV_CUSPARSE_BSR_Legacy(UCFDReal alpha, SpMat mat, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
#if defined(DEBUG)
    CheckCUDAPointer(x);
    CheckCUDAPointer(y);
#endif
    UCFDReal alpha_ = alpha, beta_ = beta;
    CUSPARSE_BSR *bsr = (CUSPARSE_BSR *)mat->data;

    CUSPARSECall(cusparseDbsrmv(
        bsr->ctx.handle, CUSPARSE_DIRECTION_ROW, CUSPARSE_OPERATION_NON_TRANSPOSE,
        bsr->mat.bn, bsr->mat.bn, bsr->bnnz, &alpha_, bsr->ctx.descr,
        ((BaseCSR *)bsr)->values, ((BaseCSR *)bsr)->rowptr,
        ((BaseCSR *)bsr)->colidx, bsr->mat.block, x, &beta_, y
    ));

    UCFDFunctionReturn(UCFD_SUCCESS);
}


extern "C" ucfd_status_t
UCFDMatCreateCUSPCSR(SpMat *mat, UCFDInt n, UCFDInt nnz, UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values)
{
    /* Check if input arrays are allocated in device side */
    CheckCUDAPointer(rowptr);
    CheckCUDAPointer(colidx);
    CheckCUDAPointer(values);

    UCFDCall(UCFDMatInit(mat));
    SpMat m = *mat;
    m->type_name = CSRCUSPARSE;
    UCFDInt alpha = 1.0;
    UCFDInt beta = 0.0;

    CUSPARSE_CSR *csr = (CUSPARSE_CSR *)calloc(1, sizeof(*csr));
    UCFDCheckNull(csr, "CSR matrix creation failed\n");

    CUSPARSECall(cusparseCreate(&csr->ctx.handle));
    CUSPARSECall(cusparseCreateCsr(
        &csr->ctx.op, n, n, nnz, rowptr, colidx, values,
        CUSPARSE_INTTYPE, CUSPARSE_INTTYPE, CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_REALTYPE
    ));

    /* Preprocess */
    CUDACall(cudaMalloc((void**)&csr->ctx.tmp, (size_t)n*sizeof(UCFDReal)));
    CUSPARSECall(cusparseCreateDnVec(&csr->ctx.vecX, n, csr->ctx.tmp, CUSPARSE_REALTYPE));
    CUSPARSECall(cusparseCreateDnVec(&csr->ctx.vecY, n, csr->ctx.tmp, CUSPARSE_REALTYPE));

    CUSPARSECall(cusparseSpMV_bufferSize(
        csr->ctx.handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, csr->ctx.op,
        csr->ctx.vecX, &beta, csr->ctx.vecY, CUSPARSE_REALTYPE, CUSPARSE_SPMV_ALG_DEFAULT,
        &csr->ctx.buffersize
    ));
    CUDACall(cudaMalloc(&csr->ctx.buffer, csr->ctx.buffersize));

#if CUSPARSE_VERSION >= 12400
    CUSPARSECall(cusparseSpMV_preprocess(
        csr->ctx.handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, csr->ctx.op,
        csr->ctx.vecX, &beta, csr->ctx.vecY, CUSPARSE_REALTYPE, CUSPARSE_SPMV_ALG_DEFAULT,
        csr->ctx.buffer
    ));
#endif

    ((BaseCSR *)csr)->n         = n;
    ((BaseCSR *)csr)->rowptr    = rowptr;
    ((BaseCSR *)csr)->colidx    = colidx;
    ((BaseCSR *)csr)->values    = values;
    m->data                     = csr;
    m->ops->spmv                = SpMV_CUSPARSE;
    m->ops->destroy             = Destroy_CUSPARSE;

    UCFDFunctionReturn(UCFD_SUCCESS);
}


extern "C" ucfd_status_t
UCFDMatCreateCUSPBSR(SpMat *mat, UCFDInt bn, UCFDInt blk, UCFDInt bnnz, UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values)
{
    /* Check if input arrays are allocated in device side */
    CheckCUDAPointer(rowptr);
    CheckCUDAPointer(colidx);
    CheckCUDAPointer(values);

    UCFDCall(UCFDMatInit(mat));
    SpMat m = *mat;
    m->type_name = BSRCUSPARSE;
    
    UCFDInt n = (UCFDInt)(bn*blk);

    CUSPARSE_BSR *bsr = (CUSPARSE_BSR *)calloc(1, sizeof(*bsr));
    UCFDCheckNull(bsr, "CSR matrix creation failed\n");
    ((BaseBSR *)bsr)->bn     = bn;
    ((BaseBSR *)bsr)->block  = blk;
    bsr->bnnz   = bnnz;

    CUSPARSECall(cusparseCreate(&bsr->ctx.handle));

#if CUSPARSE_VER_MAJOR >= 13
    UCFDInt alpha = 1.0;
    UCFDInt beta = 0.0;

    CUSPARSECall(cusparseCreateBsr(
        &bsr->ctx.op, bn, bn, bnnz, blk, blk, rowptr, colidx, values,
        CUSPARSE_INTTYPE, CUSPARSE_INTTYPE, CUSPARSE_INDEX_BASE_ZERO,
        CUSPARSE_REALTYPE, CUSPARSE_ORDER_ROW
    ));

    /* Temporal indicating of vecX and vecY */
    CUDACall(cudaMalloc((void**)&bsr->ctx.tmp, (size_t)n*sizeof(UCFDReal)));
    CUSPARSECall(cusparseCreateDnVec(&bsr->ctx.vecX, n, bsr->ctx.tmp, CUSPARSE_REALTYPE));
    CUSPARSECall(cusparseCreateDnVec(&bsr->ctx.vecY, n, bsr->ctx.tmp, CUSPARSE_REALTYPE));

    /* Preprocess */
    CUDACall(cudaMalloc((void**)&bsr->ctx.tmp, (size_t)n*sizeof(UCFDReal)));
    CUSPARSECall(cusparseCreateDnVec(&bsr->ctx.vecX, n, bsr->ctx.tmp, CUSPARSE_REALTYPE));
    CUSPARSECall(cusparseCreateDnVec(&bsr->ctx.vecY, n, bsr->ctx.tmp, CUSPARSE_REALTYPE));

    CUSPARSECall(cusparseSpMV_bufferSize(
        bsr->ctx.handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, bsr->ctx.op,
        bsr->ctx.vecX, &beta, bsr->ctx.vecY, CUSPARSE_REALTYPE, CUSPARSE_SPMV_ALG_DEFAULT,
        &bsr->ctx.buffersize
    ));
    CUDACall(cudaMalloc(&bsr->ctx.buffer, bsr->ctx.buffersize));

    CUSPARSECall(cusparseSpMV_preprocess(
        bsr->ctx.handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, bsr->ctx.op,
        bsr->ctx.vecX, &beta, bsr->ctx.vecY, CUSPARSE_REALTYPE, CUSPARSE_SPMV_ALG_DEFAULT,
        bsr->ctx.buffer
    ));
#else
    CUSPARSECall(cusparseCreateMatDescr(&bsr->ctx.descr));
#endif

    ((BaseCSR *)bsr)->n         = n;
    ((BaseCSR *)bsr)->rowptr    = rowptr;
    ((BaseCSR *)bsr)->colidx    = colidx;
    ((BaseCSR *)bsr)->values    = values;
    m->data                     = bsr;
#if CUSPARSE_VER_MAJOR >= 13
    m->ops->spmv                = SpMV_CUSPARSE;
#else
    m->ops->spmv                = SpMV_CUSPARSE_BSR_Legacy;
#endif
    m->ops->destroy             = Destroy_CUSPARSE;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

#endif
