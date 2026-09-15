#pragma once

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include <errno.h>

#include "mpisparsemat.h"


/* ------------------------------------------------------------ Internal types -- */
typedef struct {
    UCFDInt column, original;
} column_entry_t;

typedef struct {
    BASECONTEXT
    UCFDInt         *d_send_idx;
    UCFDReal        *sbuf, *ghost_x;
    UCFDReal        *d_sbuf, *d_ghost_x;
    MPI_Request     *sends, *recvs;
    int             *completed;

    int             nsend_events, peer_staging, vectorized;
    int             subgroup_width_A, subgroup_width_B;
    int             vectorized_A, vectorized_B;
    int             unique_columns;

    cudaStream_t    compute, outbound, inbound;
    cudaEvent_t     input_ready, ghost_ready, done, values_ready, *send_ready;

} CUDASpMVContext;

typedef struct {
    /* MPI context */
    CUDASpMVContext     spmvctx;

    /* Split matrix, device resident */
    BaseBSR             A, B;
    UCFDInt             *d_split_colidx;
    UCFDReal            *d_split_values;
    UCFDInt             *d_value_source;
    UCFDInt             *d_boundary_rows;
    UCFDInt             nnzb, nnzb_A, nnzb_B;
    UCFDInt             n_boundary, n_ghost;

    /* Matrix info (Host) */
    UCFDInt             block_elems;
    UCFDInt             *Arp, *Brp;
    UCFDInt             *colidx;
    UCFDInt             *value_source;
    UCFDInt             *boundary_rows;
    UCFDInt             *garray;
} CUDAMPIBSR;

/* ------------------------------------------------------------ allocation -- */

static inline size_t product(size_t a, size_t b)
{
    if (b && a > SIZE_MAX / b)
        UCFDRaiseError("allocation size overflow");
    return a * b;
}
static inline UCFDInt as_int(size_t n)
{
    if (n > INTMAX)
        
        UCFDRaiseError("32-bit BSR index / MPI count limit exceeded");
    return (UCFDInt)n;
}
static inline size_t bytes_real(size_t n)
{
    return product(n, sizeof(UCFDReal));
}
static inline size_t bytes_int(size_t n)
{
    return product(n, sizeof(UCFDInt));
}

static inline void *host_alloc(size_t count, size_t elem_size)
{
    void *p;
    product(count, elem_size);
    p = calloc(count ? count : 1, elem_size);
    if (!p)
        UCFDRaiseError("host allocation failed");
    return p;
}

#define HOST_ALLOC(type, count) ((type*)host_alloc((count), sizeof(type)))

/* Device memory allocation */
static inline void *device_alloc(size_t count, size_t elem_size)
{
    void *p = NULL;
    CUDACall(cudaMalloc(&p, product(count ? count : 1, elem_size)));
    return p;
}
static inline UCFDReal *device_alloc_real(size_t n)
{
    return (UCFDReal *)device_alloc(n, sizeof(UCFDReal));
}
static inline UCFDInt *device_alloc_int(size_t n)
{
    return (UCFDInt *)device_alloc(n, sizeof(UCFDInt));
}
static inline int grid(size_t work)
{
    size_t blocks = work / TPB + (work % TPB != 0);
    if (blocks > 65535)
        blocks = 65535; /* kernels use grid-stride loops */
    return (int)(blocks ? blocks : 1);
}
static inline UCFDInt *upload_int(const UCFDInt *src, size_t n)
{
    UCFDInt *p = device_alloc_int(n);
    if (n)
        CUDACall(cudaMemcpy(p, src, bytes_int(n), cudaMemcpyHostToDevice));
    return p;
}

static UCFDInt is_aligned16(const void *p)
{
    return (((uintptr_t)p) & 15u) == 0u;
}

static int compare_column_entry(const void *a, const void *b)
{
    const column_entry_t *x = (const column_entry_t *)a;
    const column_entry_t *y = (const column_entry_t *)b;
    if (x->column != y->column)
        return (x->column > y->column) - (x->column < y->column);
    return (x->original > y->original) - (x->original < y->original);
}

/* Full allocation ranges, not just equal bases, must be disjoint. */
static int ranges_overlap(const void *a, const void *b, size_t nbytes)
{
    uintptr_t x = (uintptr_t)a, y = (uintptr_t)b;
    return nbytes && (x <= y ? y - x < nbytes : x - y < nbytes);
}


#define DEFINE_BSR_KERNEL(NAME, BS, W, BOUNDARY, VEC)                                              \
    __global__ void NAME(int nrows, int runtime_bs, const int *__restrict__ row_list,              \
                         const int *__restrict__ rowptr, const int *__restrict__ colidx,           \
                         const double *__restrict__ a, const double *__restrict__ x, double alpha, \
                         double beta, double *__restrict__ y)                                      \
    {                                                                                              \
        const int block_size = (BS) ? (BS) : runtime_bs;                                           \
        const int lane = threadIdx.x % (W);                                                        \
        const unsigned mask = (0xffffffffu >> (32 - (W))) << ((threadIdx.x % 32) / (W) * (W));     \
        const size_t count = (size_t)nrows * block_size;                                           \
        const size_t stride = (size_t)gridDim.x * blockDim.x / (W);                                \
        for (size_t s = ((size_t)blockIdx.x * blockDim.x + threadIdx.x) / (W); s < count;          \
             s += stride)                                                                          \
        {                                                                                          \
            const int r = (BOUNDARY) ? row_list[s / block_size] : (int)(s / block_size);           \
            const int component = (int)(s % block_size);                                           \
            const int end = rowptr[r + 1];                                                         \
            double sum = 0.0;                                                                      \
            for (int k = rowptr[r]; k < end; ++k)                                                  \
            {                                                                                      \
                const double *ablk = a + ((size_t)k * block_size + component) * block_size;        \
                const double *xblk = x + (size_t)colidx[k] * block_size;                           \
                if (VEC)                                                                           \
                {                                                                                  \
                    for (int cc = 2 * lane; cc < block_size; cc += 2 * (W))                        \
                    {                                                                              \
                        const real2 av = *(const real2 *)(ablk + cc);                              \
                        const real2 xv = *(const real2 *)(xblk + cc);                              \
                        sum += av.x * xv.x + av.y * xv.y;                                          \
                    }                                                                              \
                }                                                                                  \
                else                                                                               \
                {                                                                                  \
                    for (int cc = lane; cc < block_size; cc += (W))                                \
                        sum += ablk[cc] * __ldg(xblk + cc);                                        \
                }                                                                                  \
            }                                                                                      \
            for (int offset = (W) / 2; offset; offset >>= 1)                                       \
                sum += __shfl_down_sync(mask, sum, offset, (W));                                   \
            if (!lane)                                                                             \
            {                                                                                      \
                const size_t out = (size_t)r * block_size + component;                             \
                if (BOUNDARY)                                                                      \
                    y[out] += alpha * sum;                                                         \
                else                                                                               \
                    y[out] = alpha * sum + ((beta == 0.0) ? 0.0 : beta * y[out]);                  \
            }                                                                                      \
        }                                                                                          \
    }

/* constant block sizes, scalar loads */
DEFINE_BSR_KERNEL(bsr_1_0, 1, 1, 0, 0)
DEFINE_BSR_KERNEL(bsr_1_1, 1, 1, 1, 0)
DEFINE_BSR_KERNEL(bsr_2_0, 2, 2, 0, 0)
DEFINE_BSR_KERNEL(bsr_2_1, 2, 2, 1, 0)
DEFINE_BSR_KERNEL(bsr_3_0, 3, 4, 0, 0)
DEFINE_BSR_KERNEL(bsr_3_1, 3, 4, 1, 0)
DEFINE_BSR_KERNEL(bsr_4_0, 4, 4, 0, 0)
DEFINE_BSR_KERNEL(bsr_4_1, 4, 4, 1, 0)
DEFINE_BSR_KERNEL(bsr_5_0, 5, 8, 0, 0)
DEFINE_BSR_KERNEL(bsr_5_1, 5, 8, 1, 0)
DEFINE_BSR_KERNEL(bsr_6_0, 6, 8, 0, 0)
DEFINE_BSR_KERNEL(bsr_6_1, 6, 8, 1, 0)
DEFINE_BSR_KERNEL(bsr_8_0, 8, 8, 0, 0)
DEFINE_BSR_KERNEL(bsr_8_1, 8, 8, 1, 0)
DEFINE_BSR_KERNEL(bsr_12_0, 12, 16, 0, 0)
DEFINE_BSR_KERNEL(bsr_12_1, 12, 16, 1, 0)
DEFINE_BSR_KERNEL(bsr_16_0, 16, 16, 0, 0)
DEFINE_BSR_KERNEL(bsr_16_1, 16, 16, 1, 0)
DEFINE_BSR_KERNEL(bsr_32_0, 32, 32, 0, 0)
DEFINE_BSR_KERNEL(bsr_32_1, 32, 32, 1, 0)

/* runtime block size, scalar loads */
DEFINE_BSR_KERNEL(bsr_r8_0, 0, 8, 0, 0)
DEFINE_BSR_KERNEL(bsr_r8_1, 0, 8, 1, 0)
DEFINE_BSR_KERNEL(bsr_r16_0, 0, 16, 0, 0)
DEFINE_BSR_KERNEL(bsr_r16_1, 0, 16, 1, 0)
DEFINE_BSR_KERNEL(bsr_r32_0, 0, 32, 0, 0)
DEFINE_BSR_KERNEL(bsr_r32_1, 0, 32, 1, 0)

/* constant even block sizes, 128-bit loads: W = block_size/2, one double2 per lane */
DEFINE_BSR_KERNEL(bsrv_2_0, 2, 1, 0, 1)
DEFINE_BSR_KERNEL(bsrv_2_1, 2, 1, 1, 1)
DEFINE_BSR_KERNEL(bsrv_4_0, 4, 2, 0, 1)
DEFINE_BSR_KERNEL(bsrv_4_1, 4, 2, 1, 1)
DEFINE_BSR_KERNEL(bsrv_6_0, 6, 4, 0, 1)
DEFINE_BSR_KERNEL(bsrv_6_1, 6, 4, 1, 1)
DEFINE_BSR_KERNEL(bsrv_8_0, 8, 4, 0, 1)
DEFINE_BSR_KERNEL(bsrv_8_1, 8, 4, 1, 1)
DEFINE_BSR_KERNEL(bsrv_16_0, 16, 8, 0, 1)
DEFINE_BSR_KERNEL(bsrv_16_1, 16, 8, 1, 1)
DEFINE_BSR_KERNEL(bsrv_32_0, 32, 16, 0, 1)
DEFINE_BSR_KERNEL(bsrv_32_1, 32, 16, 1, 1)

#undef DEFINE_BSR_KERNEL

/* Lanes per scalar output row. Enough to cover a block row's width without
 * splitting it across more of the warp than block_size can fill. */
static int width_for(int block_size, int vec)
{
    if (vec)
    {
        switch (block_size)
        {
        case 2:
            return 1;
        case 4:
            return 2;
        case 6:
        case 8:
            return 4;
        case 16:
            return 8;
        case 32:
            return 16;
        default:
            return 0; /* no vectorized kernel for this block_size */
        }
    }
    switch (block_size)
    {
    case 1:
        return 1;
    case 2:
        return 2;
    case 3:
    case 4:
        return 4;
    case 5:
    case 6:
    case 8:
        return 8;
    case 12:
    case 16:
        return 16;
    case 32:
        return 32;
    default:
        return block_size <= 8 ? 8 : block_size <= 16 ? 16 : 32; /* runtime fallback */
    }
}

static int scalar_count(UCFDInt nblocks, UCFDInt blocksize)
{
    return as_int(product((size_t)nblocks, (size_t)blocksize));
}

