#include <stdio.h>
#include <stdlib.h>

#include "sparsemat.h"


/**
 * CSR Matrix Format
 */
static ucfd_status_t SpMV_CSR(UCFDReal alpha, SpMat mat, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    BaseCSR *A = (BaseCSR *)mat->data;
    const UCFDInt n = A->n;
    const UCFDInt *restrict rowptr = A->rowptr, *restrict colidx = A->colidx;
    const UCFDReal *restrict values = A->values;
    UCFDInt i, j;

    OMPWrapper(j)
    for (i=0; i<n; ++i) {
        const UCFDInt st = rowptr[i];
        const UCFDInt ed = rowptr[i+1];
        UCFDReal s = 0.0;

        for (j=st; j<ed; ++j) {
            s += values[j] * x[colidx[j]];
        }
        y[i] = alpha*s + beta*y[i];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDMatCreateCSR(SpMat *mat, UCFDInt n, UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values)
{
    UCFDCall(UCFDMatInit(mat));
    SpMat m = *mat;
    m->type_name = CSR;

    BaseCSR *csr = (BaseCSR *)calloc(1, sizeof(*csr));
    UCFDCheckNull(csr, "CSR matrix creation failed\n");
    
    csr->n              = n;
    csr->rowptr         = rowptr;
    csr->colidx         = colidx;
    csr->values         = values;
    m->data             = csr;
    m->ops->spmv        = SpMV_CSR;
    m->ops->destroy     = UCFDEmptyKernel;

    UCFDFunctionReturn(UCFD_SUCCESS);
}



/** 
 * BSR Matrix format
 */
static ucfd_status_t SpMV_BSR(UCFDReal alpha, SpMat mat, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    BaseBSR *bsr = (BaseBSR *)mat->data;
    BaseCSR *A = (BaseCSR *)mat->data;
    const UCFDInt bn = bsr->bn, blk = bsr->block;
    const UCFDInt blk2 = blk*blk;
    const UCFDInt *restrict rowptr = A->rowptr, *restrict colidx = A->colidx;
    const UCFDReal *restrict values = A->values;

    UCFDInt idx, jdx, kdx, row, col;
    UCFDReal submat[blk*blk], arr[blk], xprt[blk];

    OMPWrapper(jdx, kdx, row, col, submat, arr, xprt)
    for (idx=0; idx<bn; ++idx)
    {
        // Initialize
        for (kdx=0; kdx<blk; ++kdx) arr[kdx] = 0.0;

        const UCFDInt st = rowptr[idx];
        const UCFDInt ed = rowptr[idx+1];
        
        for (jdx=st; jdx<ed; ++jdx)
        {
            const UCFDInt ldx = colidx[jdx];

            for (row=0; row<blk; ++row) {
                for (col=0; col<blk; ++col) {
                    submat[row*blk+col] = values[jdx*blk2 + row*blk + col];
                }
                xprt[row] = x[ldx*blk+row];
            }
            // blockmv
            for (row=0; row<blk; ++row) {
                UCFDReal v = 0.0;
                for (col=0; col<blk; ++col)
                    v += submat[row*blk+col] * xprt[col];
                arr[row] += v;
            }
        }

        for (kdx=0; kdx<blk; ++kdx)
            y[idx*blk+kdx] = alpha * arr[kdx] + beta*y[idx*blk+kdx];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}


ucfd_status_t UCFDMatCreateBSR(SpMat *mat, UCFDInt bn, UCFDInt blk, UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values)
{
    UCFDCall(UCFDMatInit(mat));
    SpMat m = *mat;
    m->type_name = BSR;

    BaseBSR *bsr = (BaseBSR *)calloc(1, sizeof(*bsr));
    UCFDCheckNull(bsr, "BSR matrix creation failed\n");

    bsr->bn             = bn;
    bsr->block          = blk;

    ((BaseCSR *)bsr)->n         = (UCFDInt)(bn*blk);
    ((BaseCSR *)bsr)->rowptr    = rowptr;
    ((BaseCSR *)bsr)->colidx    = colidx;
    ((BaseCSR *)bsr)->values    = values;
    m->data                     = bsr;
    m->ops->spmv                = SpMV_BSR;
    m->ops->destroy             = UCFDEmptyKernel;

    UCFDFunctionReturn(UCFD_SUCCESS);
}


#if defined(USE_MKL)
static ucfd_status_t SpMV_MKL(UCFDReal alpha, SpMat mat, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    MKLWrapper *mkl = (MKLWrapper *)mat->data;
    MKLCall(mkl_spmv(
        SPARSE_OPERATION_NON_TRANSPOSE,
        alpha, mkl->op, mkl->desc, x, beta, y
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t Destroy_MKL(SpMat mat)
{
    if (!mat) UCFDFunctionReturn(UCFD_SUCCESS);
    MKLWrapper *handle = (MKLWrapper *)mat->data;
    MKLCall(mkl_sparse_destroy(handle->op));

    UCFDFunctionReturn(UCFD_SUCCESS);    
}

ucfd_status_t UCFDMatCreateMKLBSR(SpMat *mat, UCFDInt bn, UCFDInt blk, UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values)
{
    UCFDCall(UCFDMatInit(mat));
    SpMat m = *mat;
    m->type_name = BSRMKL;

    MKLBSR *bsr = (MKLBSR *)calloc(1, sizeof(*bsr));
    UCFDCheckNull(bsr, "MKL BSR matrix creation failed\n");

    ((BaseBSR *)bsr)->bn       = bn;
    ((BaseBSR *)bsr)->block    = blk;
    bsr->handle.desc.type      = SPARSE_MATRIX_TYPE_GENERAL;
    bsr->handle.desc.mode      = 0;
    bsr->handle.desc.diag      = 0;

    MKLCall(mkl_create_bsr(
        &bsr->handle.op, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR,
        bn, bn, blk, rowptr, rowptr+1, colidx, values
    ));
    MKLCall(mkl_sparse_set_memory_hint(bsr->handle.op, MKL_MEMTYPE));

/* Optimization if `EXPECTED_SPMV_COUNT is passed */
#if defined(EXPECTED_SPMV_COUNT)
    MKLCall(mkl_sparse_set_mv_hint(
        bsr->handle.op, SPARSE_OPERATION_NON_TRANSPOSE, bsr->handle.desc, EXPECTED_SPMV_COUNT
    ));
#endif
    
    MKLCall(mkl_sparse_optimize(bsr->handle.op));

    ((BaseCSR *)bsr)->n         = (UCFDInt)(bn*blk);
    ((BaseCSR *)bsr)->rowptr    = rowptr;
    ((BaseCSR *)bsr)->colidx    = colidx;
    ((BaseCSR *)bsr)->values    = values;
    m->data                     = bsr;
    m->ops->spmv                = SpMV_MKL;
    m->ops->destroy             = Destroy_MKL;

    UCFDFunctionReturn(UCFD_SUCCESS);
}


ucfd_status_t UCFDMatCreateMKLCSR(SpMat *mat, UCFDInt n, UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values)
{
    UCFDCall(UCFDMatInit(mat));
    SpMat m = *mat;
    m->type_name = CSRMKL;

    MKLCSR *csr = (MKLCSR *)calloc(1, sizeof(*csr));
    UCFDCheckNull(csr, "MKL CSR matrix creation failed\n");

    csr->handle.desc.type      = SPARSE_MATRIX_TYPE_GENERAL;
    csr->handle.desc.mode      = 0;
    csr->handle.desc.diag      = 0;

    MKLCall(mkl_create_csr(
        &csr->handle.op, SPARSE_INDEX_BASE_ZERO, n, n, rowptr, rowptr+1, colidx, values
    ));
    MKLCall(mkl_sparse_set_memory_hint(csr->handle.op, MKL_MEMTYPE));

/* Optimization if `EXPECTED_SPMV_COUNT is passed */
#if defined(EXPECTED_SPMV_COUNT)
    MKLCall(mkl_sparse_set_mv_hint(
        bsr->handle.op, SPARSE_OPERATION_NON_TRANSPOSE, bsr->handle.desc, EXPECTED_SPMV_COUNT
    ));
#endif

    MKLCall(mkl_sparse_optimize(csr->handle.op));

    ((BaseCSR *)csr)->n         = n;
    ((BaseCSR *)csr)->rowptr    = rowptr;
    ((BaseCSR *)csr)->colidx    = colidx;
    ((BaseCSR *)csr)->values    = values;
    m->data                     = csr;
    m->ops->spmv                = SpMV_MKL;
    m->ops->destroy             = Destroy_MKL;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
#endif
