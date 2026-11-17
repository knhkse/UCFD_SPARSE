#pragma once
#include "ucfdmatimpl.h"


typedef struct {
    UCFDInt     n, nnz;
    UCFDInt     *rowptr;
    UCFDInt     *colidx;
    UCFDReal    *values;
} BaseCSR;

typedef struct {
    BaseCSR     basemat;
    UCFDInt     bn;
    UCFDInt     block;
} BaseBSR;


/* MKL matrix */
#if defined(USE_MKL)
typedef struct {
    sparse_matrix_t     op;
    struct matrix_descr desc;
} MKLWrapper;

typedef struct {
    MKLWrapper     handle;
    BaseCSR        mat;
} MKLCSR;

typedef struct {
    MKLWrapper     handle;
    BaseBSR        mat;
} MKLBSR;
#endif


/**
 * cuSPARSE matrix format
 */
#if defined(__CUDACC__)
typedef struct {
    BaseCSR         mat;
    UCFDInt         nnz;
    UCFDInt         *d_rowptr;
    UCFDInt         *d_colidx;
    cudaStream_t    stream;
} CUDACSR;

typedef struct {
    BaseBSR         mat;
    UCFDInt         nnz;
    UCFDInt         *d_rowptr;
    UCFDInt         *d_colidx;
    cudaStream_t    stream;
} CUDABSR;

typedef struct {
    cusparseHandle_t        handle;
    cusparseSpMatDescr_t    op;
    cusparseDnVecDescr_t    vecX, vecY;
    void                    *buffer;    
    size_t                  buffersize;
    UCFDReal                *tmp;
} cuSPARSEContext;

typedef struct {
    cusparseHandle_t    handle;
    cusparseMatDescr_t  descr;
} cuSPARSEContext_legacy;


typedef struct {
    CUDACSR         mat;
    cuSPARSEContext ctx;
} CUSPARSE_CSR;

typedef struct {
    BaseBSR        mat;
#if CUSPARSE_VER_MAJOR >= 13
    cuSPARSEContext ctx;
#else
    cuSPARSEContext_legacy ctx;
#endif
    UCFDInt bnnz;
} CUSPARSE_BSR;
#endif


#if defined(__cplusplus)
extern "C" {
#endif

/* Empty kernel (do nothing) */
static inline ucfd_status_t UCFDEmptyKernel(SpMat mat) {UCFDFunctionReturn(UCFD_SUCCESS);}
UCFD_INTERN ucfd_status_t UCFDMatCopyPattern(SpMat mat, UCFDInt **rp_dest, UCFDInt **ci_dest);

#if defined(__cplusplus)
}
#endif
