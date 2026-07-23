/**
 * @file        config.h
 * @brief       Header file for solver configuration
 * 
 */
#pragma once

#include <stdint.h>
#include <float.h>
#include <inttypes.h>

/**
 * Integer type designation
 */
#if defined(UCFD_INT64)
    typedef int64_t UCFDInt;
#else
    typedef int32_t UCFDInt;
#endif

/**
 * Float type designation
 */
#if defined(UCFD_FLOAT32)
    typedef float UCFDReal;
#else
    typedef double UCFDReal;
#endif

/**
 * Intel MKL configuration
 */
#if defined(USE_MKL)
    #include <mkl.h>
    #if defined(MKL_MEMTYPE_AGGRESIVE)
        #define MKL_MEMTYPE SPARSE_MEMORY_AGGRESSIVE
    #else
        #define MKL_MEMTYPE SPARSE_MEMORY_NONE
    #endif

    #if defined(UCFD_FLOAT32)
        #define mkl_create_csr mkl_sparse_s_create_csr
        #define mkl_create_bsr mkl_sparse_s_create_bsr
        #define mkl_spmv mkl_sparse_s_mv
    #else
        #define mkl_create_csr mkl_sparse_d_create_csr
        #define mkl_create_bsr mkl_sparse_d_create_bsr
        #define mkl_spmv mkl_sparse_d_mv
    #endif
#endif

/**
 * CUDA configuration
 * Import cuBLAS automatically
 * Currently, Krylov subspace methods for CUDA computation
 * require cuBLAS functions
 */
#if defined(USE_CUDA)
    #include <cuda_runtime.h>
    #include <cublas_v2.h>

    #if !defined(TPB)
        #define TPB 128         // Default Threads-per-block size
    #endif
#endif

/**
 * cuSPARSE configuration
 */
#if defined(USE_CUSPARSE)
    #include <cusparse.h>

    #if defined(UCFD_INT64)
        #define CUSPARSE_INTTYPE CUSPARSE_INDEX_64I
    #else
        #define CUSPARSE_INTTYPE CUSPARSE_INDEX_32I
    #endif

    #if defined(UCFD_FLOAT32)
        #define CUSPARSE_REALTYPE CUDA_R_32F
    #else
        #define CUSPARSE_REALTYPE CUDA_R_64F
    #endif
#endif
