#pragma once

#include <stdio.h>
#include <stdlib.h>

/**
 * Vendor headers (mkl.h, cuda_runtime.h, cublas_v2.h, cusparse.h) are owned
 * by config.h. This file must be reachable only through ucfd_types.h (which
 * includes config.h first), never included on its own or before config.h -
 * otherwise the *Call macros below reference vendor status types that don't
 * exist yet.
 */

#define UCFDFunctionReturn(...) return __VA_ARGS__
#define UCFD_EXTERN extern __attribute__((visibility("default")))
#define UCFD_INTERN extern __attribute__((visibility("hidden")))
#define MAYBE_UNUSED __attribute__((unused))

#if defined(USE_OMP)
    #define pragma_indvars(x) _Pragma(#x)
    #define OMPWrapper(...) pragma_indvars(omp parallel for private(__VA_ARGS__))
    #define OMPSumReduction(...) pragma_indvars(omp parallel for reduction(+:__VA_ARGS__))
    #define OMPFOR pragma_indvars(omp parallel for)
    #define OMPScheduleStaticSumReduction(...) pragma_indvars(omp parallel for schedule(static) reduction(+:__VA_ARGS__))
#else
    #define OMPWrapper(...)
    #define OMPSumReduction(...)
    #define OMPFOR
    #define OMPScheduleStaticSumReduction(...)
#endif

#define UCFDWarning(msg) fprintf(stderr, "%s", msg);

#define UCFDCheckNull(obj, msg)                                             \
    do {                                                                    \
        if (!(obj)) {                                                       \
            fprintf(stderr, msg);                                           \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    } while (0)

#define UCFDMatch(obj1, obj2, ...)                                          \
    do {                                                                    \
        if (obj1 != obj2) fprintf(stderr, __VA_ARGS__);                     \
    } while (0)


/**
 * Library error-checking macros
 */
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

#define CheckCUDAPointer(ptr)                                               \
    do {                                                                    \
        UCFDCheckNull(ptr, "Null pointer is given\n");                      \
        struct cudaPointerAttributes attr;                                  \
        CUDACall(cudaPointerGetAttributes(&attr, ptr));                     \
        if (attr.type != cudaMemoryTypeDevice) {                            \
            fprintf(stderr, "%s:%d : Device pointer must be passed\n",      \
                    __FILE__, __LINE__);                                    \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    } while (0)
