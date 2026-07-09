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


