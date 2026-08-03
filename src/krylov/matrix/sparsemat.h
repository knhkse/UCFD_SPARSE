#pragma once
#include "ucfdmatimpl.h"


typedef struct {
    UCFDInt     n;
    UCFDInt     *rowptr;
    UCFDInt     *colidx;
    UCFDReal    *values;
} BaseCSR;

typedef struct {
    BaseCSR    basemat;
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
    BaseCSR        mat;
    MKLWrapper     handle;
} MKLCSR;

typedef struct {
    BaseBSR        mat;
    MKLWrapper     handle;
} MKLBSR;
#endif

/**
 * cuSPARSE matrix format
 */
#if defined(USE_CUSPARSE)
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
    BaseCSR        mat;
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


/* Empty kernel (do nothing) */
static inline ucfd_status_t UCFDEmptyKernel(SpMat mat) {UCFDFunctionReturn(UCFD_SUCCESS);}
