#pragma once
#include "ucfdmatimpl.h"


/**
 * BSR matrix
 */
typedef struct {
    UCFDInt bn;
    UCFDInt block;
} SpMat_BSR;



/**
 * CSR matrix
 */
typedef struct {
    char dummy;     /* Dummy component */
} SpMat_CSR;










/**
 * MKL matrix format
 */
#if defined(USE_MKL)
typedef struct {
    sparse_matrix_t op;
    struct matrix_descr desc;
} MKLWrapper;

typedef struct {
    MKLWrapper handle;
} SpMat_MKLCSR;

typedef struct {
    MKLWrapper handle;
    UCFDInt bn;
    UCFDInt block;
} SpMat_MKLBSR;
#endif

/**
 * cuSPARSE matrix format
 */
#if defined(USE_CUSPARSE)
typedef struct {
    cusparseHandle_t handle;
    cusparseSpMatDescr_t op;
    cusparseDnVecDescr_t vecX, vecY;
    void *buffer;
    size_t buffersize;
    UCFDReal *tmp;
} cuSPARSEContext;

typedef struct {
    cusparseHandle_t handle;
    cusparseMatDescr_t descr;
} cuSPARSEContext_legacy;


typedef struct {
    cuSPARSEContext ctx;
} SpMat_CUSPARSECSR;

typedef struct {
#if CUSPARSE_VER_MAJOR >= 13
    cuSPARSEContext ctx;
#else
    cuSPARSEContext_legacy ctx;
#endif
    UCFDInt bn;
    UCFDInt block;
    UCFDInt bnnz;
} SpMat_CUSPARSEBSR;
#endif

/* Empty kernel (do nothing) */
static inline ucfd_status_t UCFDEmptyKernel(SpMat mat) {UCFDFunctionReturn(UCFD_SUCCESS);}
