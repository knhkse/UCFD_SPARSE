#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "ucfd_types.h"

#if defined(USE_MKL)
    #include <mkl.h>
#endif

#if defined(USE_CUDA)
    #include <cuda_runtime.h>
#endif

#if defined(USE_CUBLAS)
    #include <cublas_v2.h>
#endif

#if defined(USE_CUSPARSE)
    #include <cusparse.h>
#endif


#define UCFDFunctionReturn(...) return __VA_ARGS__
#define UCFD_EXTERN extern __attribute__((visibility("default")))
#define UCFD_INTERN extern __attribute__((visibility("hidden")))

#if defined(USE_OMP)
    #define pragma_indvars(x) _Pragma(#x)
    #define OMPWrapper(...) pragma_indvars(omp parallel for private(__VA_ARGS__))
    #define OMPSumReduction(...) pragma_indvars(omp parallel for reduction(+:__VA_ARGS__))
#else
    #define OMPWrapper(...)
    #define OMPSumReduction(...)
#endif


#define UCFDNullCheck(obj, msg)                                             \
    do {                                                                    \
        if (!(obj)) fprintf(stderr, msg);                                   \
    } while (0)

#define UCFDCompare(obj1, obj2, ...)                                        \
    do {                                                                    \
        if (obj1 != obj2) fprintf(stderr, __VA_ARGS__);                     \
    } while (0)

#define UCFDError(msg)                                                      \
    do {                                                                    \
        fprintf(stderr, msg);                                               \
        exit(EXIT_FAILURE);                                                 \
    } while (0)


#define UCFDCall(call)                                                      \
    do {                                                                    \
        ucfd_status_t st_ = (call);                                         \
        if (st_ != UCFD_SUCCESS) {                                          \
            fprintf(stderr, "[UCFD] %s:%d  status %d\n",                    \
                    __FILE__, __LINE__, (int)st_);                          \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    } while (0)

#define MKLCall(call)                                                       \
    do {                                                                    \
        sparse_status_t st_ = (call);                                       \
        if (st_ != SPARSE_STATUS_SUCCESS) {                                 \
            fprintf(stderr, "[MKL]   %s:%d  status %d\n",                   \
                    __FILE__, __LINE__, (int)st_);                          \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    } while (0)

#define CUDACall(call)                                                      \
    do {                                                                    \
        cudaError_t err_ = (call);                                          \
        if (err_ != cudaSuccess) {                                          \
            fprintf(stderr, "[CUDA]   %s:%d  %s\n",                         \
                    __FILE__, __LINE__, cudaGetErrorString(err_));          \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    } while (0)

#define CUBLASCall(call)                                                    \
    do {                                                                    \
        cublasStatus_t st_ = (call);                                        \
        if (st_ != CUBLAS_STATUS_SUCCESS) {                                 \
            fprintf(stderr, "[cuBLAS] %s:%d  status %d\n",                  \
                    __FILE__, __LINE__, (int)st_);                          \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    } while (0)

#define CUSPARSECall(call)                                                  \
    do {                                                                    \
        cusparseStatus_t st_ = (call);                                      \
        if (st_ != CUSPARSE_STATUS_SUCCESS) {                               \
            fprintf(stderr, "[cuSPARSE] %s:%d  %s\n",                       \
                    __FILE__, __LINE__, cusparseGetErrorString(st_));       \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    } while (0)

