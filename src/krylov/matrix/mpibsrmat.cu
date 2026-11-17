#include "mpihelper.h"
#include "mpisparsemat.cuh"


/* --------------------------------------------------------- steady kernels */
/* The halo moves whole blocks: one index selects block_size consecutive components */
__global__ void
gather_kernel(size_t n, UCFDInt block_size, const UCFDInt *__restrict__ indices,
              const UCFDReal *__restrict__ x, UCFDReal *__restrict__ packed)
{
    for (size_t s = (size_t)blockIdx.x * blockDim.x + threadIdx.x; s < n;
         s += (size_t)gridDim.x * blockDim.x)
        packed[s] = x[(size_t)indices[s / block_size] * block_size + s % block_size];
}

__global__ void
scale_kernel(size_t n, UCFDReal beta, UCFDReal *__restrict__ y)
{
    for (size_t s = (size_t)blockIdx.x * blockDim.x + threadIdx.x; s < n;
         s += (size_t)gridDim.x * blockDim.x)
        y[s] = (beta == 0.0) ? 0.0 : beta * y[s]; /* never read old y at beta=0 */
}

__global__ void
permute_values_kernel(size_t n, size_t block_elems,
                      const UCFDInt *__restrict__ value_source,
                      const UCFDReal *__restrict__ src, UCFDReal *__restrict__ dst)
{
    for (size_t s = (size_t)blockIdx.x * blockDim.x + threadIdx.x; s < n;
         s += (size_t)gridDim.x * blockDim.x)
        dst[s] = __ldg(&src[(size_t)value_source[s / block_elems] * block_elems + s % block_elems]);
}


/* ------------------------------------------------------------- matvec kernels */
static void bsr_mult_custom(CUDASpMVContext *ctx, int boundary_only, int vectorized,
                            int subgroup_width, UCFDInt nblock_rows, const UCFDInt *row_list,
                            BaseBSR *matrix, const UCFDReal *x, UCFDReal alpha, UCFDReal beta, UCFDReal *y)
{
    const UCFDInt block_size = ctx->blocksize;
    const int blocks =
        grid(product(product((size_t)nblock_rows, (size_t)block_size), (size_t)subgroup_width));
#define ARGS                                                                                       \
    nblock_rows, block_size, row_list, matrix->basemat.rowptr, matrix->basemat.colidx,             \
    matrix->basemat.values, x, alpha, beta, y
#define PICK(plain, bnd)                                                                           \
    do                                                                                             \
    {                                                                                              \
        if (boundary_only)                                                                         \
            bnd<<<blocks, TPB, 0, ctx->compute>>>(ARGS);                                           \
        else                                                                                       \
            plain<<<blocks, TPB, 0, ctx->compute>>>(ARGS);                                         \
    } while (0)

    if (vectorized)
    {
        switch (block_size)
        {
        case 2:
            PICK(bsrv_2_0, bsrv_2_1);
            break;
        case 4:
            PICK(bsrv_4_0, bsrv_4_1);
            break;
        case 6:
            PICK(bsrv_6_0, bsrv_6_1);
            break;
        case 8:
            PICK(bsrv_8_0, bsrv_8_1);
            break;
        case 16:
            PICK(bsrv_16_0, bsrv_16_1);
            break;
        case 32:
            PICK(bsrv_32_0, bsrv_32_1);
            break;
        default:
            UCFDRaiseError("no vectorized kernel for this block size");
        }
    }
    else
    {
        switch (block_size)
        {
        case 1:
            PICK(bsr_1_0, bsr_1_1);
            break;
        case 2:
            PICK(bsr_2_0, bsr_2_1);
            break;
        case 3:
            PICK(bsr_3_0, bsr_3_1);
            break;
        case 4:
            PICK(bsr_4_0, bsr_4_1);
            break;
        case 5:
            PICK(bsr_5_0, bsr_5_1);
            break;
        case 6:
            PICK(bsr_6_0, bsr_6_1);
            break;
        case 8:
            PICK(bsr_8_0, bsr_8_1);
            break;
        case 12:
            PICK(bsr_12_0, bsr_12_1);
            break;
        case 16:
            PICK(bsr_16_0, bsr_16_1);
            break;
        case 32:
            PICK(bsr_32_0, bsr_32_1);
            break;
        default:
            if (subgroup_width == 8)
                PICK(bsr_r8_0, bsr_r8_1);
            else if (subgroup_width == 16)
                PICK(bsr_r16_0, bsr_r16_1);
            else
                PICK(bsr_r32_0, bsr_r32_1);
            break;
        }
    }
#undef PICK
#undef ARGS
    CUDACall(cudaGetLastError());
}

static void bsr_mult(CUDAMPIBSR *mat, UCFDReal alpha, const UCFDReal *d_x_local, UCFDReal beta, UCFDReal *d_y)
{
    CUDASpMVContext *ctx = &mat->spmvctx;

    const size_t local = product((size_t)mat->A.bn, (size_t)ctx->blocksize);
    int vec;
    if (!local) return;
    if (alpha == 0.0 || mat->nnzb_A == 0)
    {
        /* beta must still be applied exactly once, including on empty nblock_rows. */
        if (beta != 1.0)
        {
            scale_kernel<<<grid(local), TPB, 0, ctx->compute>>>(local, beta, d_y);
            CUDACall(cudaGetLastError());
        }
        return;
    }
    vec = ctx->vectorized_A && is_aligned16(d_x_local);
    bsr_mult_custom(
        ctx, 0, vec, width_for(ctx->blocksize, vec), mat->A.bn, NULL,
        &mat->A, d_x_local, alpha, beta, d_y
    );
}

static ucfd_status_t SpMV_CUDAMPIBSR(UCFDReal alpha, SpMat mat, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    /**
     * Notes
     * x and y are device arrays;
     */
    CUDAMPIBSR *bsr = (CUDAMPIBSR *)mat->data;
    CUDASpMVContext *ctx = &bsr->spmvctx;

    const UCFDInt block_size = ctx->blocksize;
    const size_t send_scalars = product((size_t)ctx->send_total, (size_t)block_size);
    int next_send = 0, received = 0;

#if defined(DEBUG)
    if (bn && (!x || !y || ranges_overlap(x, y, bytes_real((size_t)bn * block_size))))
        UCFDRaiseError("SpMV needs distinct device x and y");
#endif

    if (!ctx->nsend && !ctx->nrecv) {
        /* no peers : no host waits at all */
        bsr_mult(bsr, alpha, x, beta, y);
        CUDACall(cudaEventRecord(ctx->done, ctx->compute));
        UCFDFunctionReturn(UCFD_SUCCESS);
    }

    /* The previous upload must be finished before MPI overwrites ghost_x, and
     * it trails the previous compute, so this also orders old ghost readers
     * ahead of the new writes. */
    CUDACall(cudaEventSynchronize(ctx->ghost_ready));
    if (ctx->nrecv) MPICall(MPI_Startall(ctx->nrecv, ctx->recvs));
    CUDACall(cudaEventRecord(ctx->input_ready, ctx->compute));
    CUDACall(cudaStreamWaitEvent(ctx->outbound, ctx->input_ready, 0));
    CUDACall(cudaStreamWaitEvent(ctx->inbound, ctx->input_ready, 0));
    
    if (send_scalars)
    {
        /* Gather on the GPU: only the blocks peers actually asked for cross
         * PCIe, never the whole local vector. */
        gather_kernel<<<grid(send_scalars), TPB, 0, ctx->outbound>>>(
            send_scalars, block_size, ctx->d_send_idx, x, ctx->d_sbuf
        );
        CUDACall(cudaGetLastError());
        if (ctx->peer_staging) {
            int j;
            for (j=0; j<ctx->nsend; ++j)
            {
                const size_t off = (size_t)ctx->send_off[j] * block_size;
                CUDACall(cudaMemcpyAsync(
                    ctx->sbuf + off, ctx->d_sbuf + off,
                    bytes_real((size_t)ctx->send_count[j] * block_size),
                    cudaMemcpyDeviceToHost, ctx->outbound
                ));
                CUDACall(cudaEventRecord(ctx->send_ready[j], ctx->outbound));
            }
        } else {
            CUDACall(cudaMemcpyAsync(ctx->sbuf, ctx->d_sbuf, bytes_real(send_scalars),
                                     cudaMemcpyDeviceToHost, ctx->outbound));
            CUDACall(cudaEventRecord(ctx->send_ready[0], ctx->outbound));
        }
    }

    // Hides halo
    bsr_mult(bsr, alpha, x, beta, y);

    while (next_send < ctx->nsend || received < ctx->nrecv)
    {
        /* All D2H copies share one stream, so finished peers form a prefix:
         * poll only the next event instead of rescanning every peer. */
        while (next_send < ctx->nsend)
        {
            cudaError_t state = cudaEventQuery(ctx->send_ready[ctx->peer_staging ? next_send : 0]);
            if (state == cudaErrorNotReady)
                break;
            CUDACall(state);
            if (ctx->peer_staging) MPICall(MPI_Start(&ctx->sends[next_send++]));
            else {
                MPICall(MPI_Startall(ctx->nsend, ctx->sends));
                next_send = ctx->nsend;
            }
        }
        if (received < ctx->nrecv)
        {
            int count = 0, q;
            MPICall(
                MPI_Testsome(ctx->nrecv, ctx->recvs, &count, ctx->completed, MPI_STATUSES_IGNORE
            ));
            if (count != MPI_UNDEFINED)
                for (q = 0; q < count; ++q)
                {
                    const int j = ctx->completed[q];
                    const size_t off = (size_t)ctx->recv_off[j] * block_size;
                    if (ctx->peer_staging)
                        CUDACall(cudaMemcpyAsync(
                            ctx->d_ghost_x + off, ctx->ghost_x + off,
                            bytes_real((size_t)ctx->recv_count[j] * block_size),
                            cudaMemcpyHostToDevice, ctx->inbound
                        ));
                    ++received;
                }
        }
        else if (next_send < ctx->nsend)
        {
            /* A rank that only sends still has to enter MPI so its transport
             * makes progress while it waits on CUDA. */
            int flag;
            MPICall(MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, ctx->comm, &flag, MPI_STATUS_IGNORE));
        }
    }
    if (!ctx->peer_staging && ctx->nrecv)
        CUDACall(cudaMemcpyAsync(
            ctx->d_ghost_x, ctx->ghost_x,
            bytes_real(product((size_t)bsr->n_ghost, (size_t)block_size)),
            cudaMemcpyHostToDevice, ctx->inbound
        ));
    CUDACall(cudaEventRecord(ctx->ghost_ready, ctx->inbound));

    if (alpha != 0.0 && bsr->n_boundary)
    {
        /* Both writers of the boundary rows of y sit on one stream, so the
         * accumulate cannot race the owned product. */
        CUDACall(cudaStreamWaitEvent(ctx->compute, ctx->ghost_ready, 0));
        bsr_mult_custom(ctx, 1, ctx->vectorized_B, ctx->subgroup_width_B,
                        bsr->n_boundary, bsr->d_boundary_rows, &bsr->B,
                        ctx->d_ghost_x, alpha, 1.0, y);
    }
    CUDACall(cudaEventRecord(ctx->done, ctx->compute));
    /* Sends matter only for buffer reuse, not for the result, so waiting here
     * lets the boundary_block_rows kernel overlap the outstanding sends. */
    if (ctx->nsend) MPICall(MPI_Waitall(ctx->nsend, ctx->sends, MPI_STATUSES_IGNORE));

    UCFDSpMVSynchronize(ctx);

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t inline UCFDSpMVSynchronize(CUDASpMVContext *ctx)
{
    CUDACall(cudaEventSynchronize(ctx->done));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

/* =========================================================================
 * Host-side A/B split
 *
 * Structure only, in block space: the same two-pass analysis split_csr()
 * performs, with a block index everywhere a scalar index used to be. Values
 * are not read, not copied and not even mapped -- only `value_source`, which records
 * for each split block which original block it came from.
 *
 * Layout of the split, matching the device arrays exactly:
 *   blocks [0, nnzb_A)      A, owned columns, remapped to local 0-based
 *   blocks [nnzb_A, nnzb)   B, ghost columns, remapped to a garray position
 * ========================================================================= */
static ucfd_status_t UCFDSplitBSR(UCFDInt bn, UCFDInt blk,
                                  const UCFDInt *rp, const UCFDInt *ci,
                                  UCFDInt cstart, UCFDInt cend,
                                  CUDAMPIBSR *mat)
{
    UCFDInt *remote;
    column_entry_t *order;
    UCFDInt nremote = 0, max_row = 0;
    mat->spmvctx.unique_columns = 1;

#if defined(DEBUG)
    if (!rp || rp[0] != 0)
        UCFDRaiseError("Host rp must start at zero, including empty ranks");
#endif

    for (UCFDInt ib = 0; ib < bn; ++ib)
    {
        UCFDInt length;
        if (rp[ib] < 0 || rp[ib + 1] < rp[ib])
            UCFDRaiseError("rp must be nonnegative and nondecreasing");
        length = rp[ib + 1] - rp[ib];
        if (length > max_row)
            max_row = length;
    }
    mat->nnzb = rp[bn];
    bytes_real(product((size_t)mat->nnzb, blk*blk));
    mat->Arp = HOST_ALLOC(UCFDInt, (size_t)bn + 1);
    mat->Brp = HOST_ALLOC(UCFDInt, (size_t)bn + 1);
    remote = HOST_ALLOC(UCFDInt, mat->nnzb);

    for (UCFDInt ib = 0; ib < bn; ++ib)
    {
        UCFDInt na = 0;
        for (UCFDInt kb = rp[ib]; kb < rp[ib + 1]; ++kb)
        {
            const UCFDInt g = ci[kb];
            if (g >= cstart && g < cend)
                ++na;
            else
                remote[nremote++] = g;
        }
        mat->Arp[ib + 1] = na;
        mat->Brp[ib + 1] = rp[ib + 1] - rp[ib] - na;
    }

    /* Sorted unique global block columns define compact ghost block slots. */
    qsort(remote, nremote, sizeof(*remote), cmp_idx);
    for (UCFDInt k = 0; k < nremote; ++k)
        if (!k || remote[k] != remote[k - 1])
            remote[mat->n_ghost++] = remote[k];
    mat->garray = HOST_ALLOC(UCFDInt, mat->n_ghost);
    memcpy(mat->garray, remote, bytes_int(mat->n_ghost));
    free(remote);

    for (UCFDInt ib = 0; ib < bn; ++ib)
    {
        mat->Arp[ib + 1] += mat->Arp[ib];
        mat->Brp[ib + 1] += mat->Brp[ib];
    }
    mat->nnzb_A = mat->Arp[bn];
    mat->nnzb_B = mat->Brp[bn];
    mat->colidx = HOST_ALLOC(UCFDInt, mat->nnzb);
    mat->value_source = HOST_ALLOC(UCFDInt, mat->nnzb);
    order = HOST_ALLOC(column_entry_t, max_row);
    for (UCFDInt ib = 0; ib < bn; ++ib)
    {
        const UCFDInt length = rp[ib + 1] - rp[ib];
        UCFDInt ap = mat->Arp[ib], bp = mat->nnzb_A + mat->Brp[ib];
        int sorted = 1;
        for (UCFDInt j = 0; j < length; ++j)
        {
            order[j].column = ci[rp[ib] + j];
            order[j].original = rp[ib] + j;
            if (j && order[j - 1].column > order[j].column)
                sorted = 0;
        }
        if (!sorted)
            qsort(order, length, sizeof(*order), compare_column_entry);
        for (UCFDInt j = 0; j < length; ++j)
        {
            const UCFDInt g = order[j].column, kb = order[j].original;
            UCFDInt dest;
            if (j && order[j - 1].column == g)
                mat->spmvctx.unique_columns = 0;
            if (g >= cstart && g < cend)
            {
                dest = ap++;
                mat->colidx[dest] = g - cstart;
            }
            else
            {
                dest = bp++;
                mat->colidx[dest] = bsearch_idx(mat->garray, mat->n_ghost, g);
                if (mat->colidx[dest] < 0)
                    UCFDRaiseError("missing ghost block in garray");
            }
            mat->value_source[dest] = kb;
        }
    }
    free(order);
    for (UCFDInt ib = 0; ib < bn; ++ib)
        mat->n_boundary += mat->Brp[ib + 1] != mat->Brp[ib];
    mat->boundary_rows = HOST_ALLOC(UCFDInt, mat->n_boundary);
    for (UCFDInt ib = 0, q = 0; ib < bn; ++ib)
        if (mat->Brp[ib + 1] != mat->Brp[ib])
            mat->boundary_rows[q++] = ib;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static void select_kernels(UCFDInt blocksize, CUDAMPIBSR *mat)
{
    CUDASpMVContext *ctx = &mat->spmvctx;
    const int wide = ctx->vectorized && width_for(blocksize, 1) > 0;
    ctx->vectorized_A = wide && is_aligned16(mat->A.basemat.values);
    ctx->vectorized_B = wide && is_aligned16(mat->B.basemat.values) && \
                        is_aligned16(ctx->d_ghost_x);
    ctx->subgroup_width_A = width_for(blocksize, ctx->vectorized_A);
    ctx->subgroup_width_B = width_for(blocksize, ctx->vectorized_B);
}

static void apply_values(CUDAMPIBSR *mat, const UCFDReal *d_val)
{
    const size_t n = product((size_t)mat->nnzb, (size_t)mat->block_elems);
    
    if (n)
    {
    #if defined(DEBUG)
        if (!d_val || ranges_overlap(d_val, mat->d_split_values, bytes_double(n)))
            UCFDRaiseError("value source must be a separate device array");
    #endif
        permute_values_kernel<<<grid(n), TPB, 0, mat->spmvctx.compute>>>(
            n, mat->block_elems, mat->d_value_source, d_val, mat->d_split_values
        );
        CUDACall(cudaGetLastError());
    }
    CUDACall(cudaEventRecord(mat->spmvctx.done, mat->spmvctx.compute));
}

/* build_schedule */
static ucfd_status_t UCFDSetCUDAMPIBSRContext(CUDASpMVContext *ctx,
                                              const UCFDInt *block_range,
                                              const UCFDInt *garray,
                                              UCFDInt n_ghost,
                                              UCFDInt blocksize)
{
    const int size = ctx->size;
    int *req_to = HOST_ALLOC(int, size), *req_from = HOST_ALLOC(int, size);
    int *req_local = HOST_ALLOC(int, n_ghost);
    int p, j, k, off, r;
    MPI_Request *rq;

    for (k = 0; k < n_ghost; ++k)
    {
        p = owner(block_range, size, garray[k]);
        req_local[k] = garray[k] - block_range[p];
        req_to[p]++;
    }
    ctx->nrecv = 0;
    for (p = 0; p < size; ++p)
        if (req_to[p])
            ctx->nrecv++;
    ctx->recv_nei = HOST_ALLOC(int, ctx->nrecv);
    ctx->recv_count = HOST_ALLOC(int, ctx->nrecv);
    ctx->recv_off = HOST_ALLOC(int, ctx->nrecv);
    for (p = 0, j = 0, off = 0; p < size; ++p)
        if (req_to[p])
        {
            ctx->recv_nei[j] = p;
            ctx->recv_count[j] = req_to[p];
            ctx->recv_off[j] = off;
            off += req_to[p];
            ++j;
        }
    MPICall(MPI_Alltoall(req_to, 1, MPI_INT, req_from, 1, MPI_INT, ctx->comm));

    ctx->nsend = 0;
    for (p = 0; p < size; ++p)
        if (req_from[p])
            ctx->nsend++;
    ctx->send_nei = HOST_ALLOC(int, ctx->nsend);
    ctx->send_count = HOST_ALLOC(int, ctx->nsend);
    ctx->send_off = HOST_ALLOC(int, ctx->nsend);
    ctx->send_total = 0;
    for (p = 0, j = 0; p < size; ++p)
        if (req_from[p])
        {
            ctx->send_nei[j] = p;
            ctx->send_count[j] = req_from[p];
            ctx->send_off[j] = ctx->send_total;
            ctx->send_total = as_int((size_t)ctx->send_total + req_from[p]);
            ++j;
        }
    ctx->send_idx = HOST_ALLOC(int, ctx->send_total);

    rq = HOST_ALLOC(MPI_Request, as_int((size_t)ctx->nrecv + ctx->nsend));
    r = 0;
    for (j = 0; j < ctx->nsend; ++j)
        MPICall(MPI_Irecv(ctx->send_idx + ctx->send_off[j], ctx->send_count[j], MPI_INT,
                         ctx->send_nei[j], ctx->tag, ctx->comm, &rq[r++]));
    for (j = 0; j < ctx->nrecv; ++j)
        MPICall(MPI_Isend(req_local + ctx->recv_off[j], ctx->recv_count[j], MPI_INT,
                         ctx->recv_nei[j], ctx->tag, ctx->comm, &rq[r++]));
    if (r) MPICall(MPI_Waitall(r, rq, MPI_STATUSES_IGNORE));
    free(rq);
    free(req_local);
    free(req_to);
    free(req_from);

    /* Every halo message is block_size doubles per block; check the MPI count fits. */
    for (j = 0; j < ctx->nsend; ++j)
        (void)as_int(product((size_t)ctx->send_count[j], (size_t)blocksize));
    for (j = 0; j < ctx->nrecv; ++j)
        (void)as_int(product((size_t)ctx->recv_count[j], (size_t)blocksize));

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static void bsr_spmv_synchronize(CUDASpMVContext *ctx)
{
    CUDACall(cudaEventSynchronize(ctx->done));
}

static ucfd_status_t UCFDMPIMatDestroy(SpMat mat)
{
    if (!mat) UCFDFunctionReturn(UCFD_SUCCESS);
    CUDAMPIBSR *bsr = (CUDAMPIBSR *)mat->data;
    CUDASpMVContext *ctx = &bsr->spmvctx;

    /* Context destroy */
    int j;
    CUDACall(cudaStreamSynchronize(ctx->compute));
    CUDACall(cudaStreamSynchronize(ctx->outbound));
    CUDACall(cudaStreamSynchronize(ctx->inbound));
    for (j=0; j<ctx->nrecv; ++j)
        MPICall(MPI_Request_free(&ctx->recvs[j]));
    for (j=0; j<ctx->nsend; ++j)
        MPICall(MPI_Request_free(&ctx->sends[j]));
    for (j=0; j<ctx->nsend_events; ++j)
        CUDACall(cudaEventDestroy(ctx->send_ready[j]));

    CUDACall(cudaEventDestroy(ctx->input_ready));
    CUDACall(cudaEventDestroy(ctx->ghost_ready));
    CUDACall(cudaEventDestroy(ctx->done));
    CUDACall(cudaEventDestroy(ctx->values_ready));
    CUDACall(cudaStreamDestroy(ctx->compute));
    CUDACall(cudaStreamDestroy(ctx->outbound));
    CUDACall(cudaStreamDestroy(ctx->inbound));
    CUDACall(cudaFreeHost(ctx->sbuf));
    CUDACall(cudaFreeHost(ctx->ghost_x));
    CUDACall(cudaFree(ctx->d_send_idx));
    CUDACall(cudaFree(ctx->d_sbuf));
    CUDACall(cudaFree(ctx->d_ghost_x));
    free(ctx->sends);
    free(ctx->recvs);
    free(ctx->completed);
    free(ctx->send_ready);
    free(ctx->recv_nei);
    free(ctx->recv_count);
    free(ctx->recv_off);
    free(ctx->send_nei);
    free(ctx->send_count);
    free(ctx->send_off);
    free(ctx->send_idx);
    ctx->comm = MPI_COMM_NULL;

    /* BSR matrix destroy */
    free(bsr->Arp);
    free(bsr->colidx);
    CUDACall(cudaFree(bsr->A.basemat.rowptr));
    CUDACall(cudaFree(bsr->B.basemat.rowptr));
    CUDACall(cudaFree(bsr->d_split_colidx));
    CUDACall(cudaFree(bsr->d_split_values));
    CUDACall(cudaFree(bsr->d_value_source));
    CUDACall(cudaFree(bsr->d_boundary_rows));

    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDMatUpdateDeviceValues(SpMat mat, UCFDReal *d_values)
{
    CUDAMPIBSR *bsr = (CUDAMPIBSR *)mat->data;
    apply_values(bsr, d_values);

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t UCFDCopyComputeStream(SpMat mat, cudaStream_t *stream)
{
    CUDAMPIBSR *bsr = (CUDAMPIBSR *)mat->data;

    *stream = bsr->spmvctx.compute;
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t UCFDBSRCopyPattern(SpMat mat, UCFDInt **rp_dest, UCFDInt **ci_dest)
{
    CUDAMPIBSR *bsr = (CUDAMPIBSR *)mat->data;
    *rp_dest = bsr->Arp;
    *ci_dest = bsr->colidx;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t UCFDBSRCopyValues(SpMat mat, UCFDReal *__restrict__ values)
{
    CUDAMPIBSR *bsr = (CUDAMPIBSR *)mat->data;
    const UCFDInt val_count = product((size_t)bsr->nnzb_A, bsr->block_elems);
    CUDACall(cudaMemcpyAsync(
        values, bsr->d_split_values, val_count*sizeof(UCFDReal),
        cudaMemcpyDeviceToDevice, bsr->spmvctx.compute
    ));

    UCFDFunctionReturn(UCFD_SUCCESS);
}

/* Collective creation. rp/ci are HOST arrays; d_va is DEVICE storage in
 * original BSR order. va_stream identifies its producer, including stream 0.
 * Explicit event ordering is required because internal streams are nonblocking.
 * Creation waits for the initial permutation, so all input arrays may be freed
 * after return. The context owns everything needed for later SpMVs/updates.
 * block_range may be NULL (derive it), or a checked common ownership prefix.
 * Select the rank's GPU before calling; destroy before MPI_Finalize().
 */
extern "C" ucfd_status_t
UCFDMatCreateCUDAMPIBSR(SpMat *mat, Ctx ctx,
                        UCFDInt bn, UCFDInt blk,
                        UCFDInt *rp, UCFDInt *ci, UCFDReal *d_va,
                        UCFDInt peer_staging, UCFDInt vectorize)
{
    int j, low_priority, high_priority;

    UCFDCall(UCFDMatInit(mat));
    SpMat m = *mat;
    m->type_name = BSRCUDAMPI;

    CUDAMPIBSR *bsr = (CUDAMPIBSR *)calloc(1, sizeof(*bsr));
    UCFDCheckNull(bsr, "CUDAMPIBSR matrix creation failed\n");

    CUDASpMVContext *c = &bsr->spmvctx;

    ContextNextTag(ctx, &c->tag);
    c->comm = ctx->comm;

    int size=0, rank=0;
    MPICall(MPI_Comm_set_errhandler(ctx->comm, MPI_ERRORS_RETURN));
    MPICall(MPI_Comm_rank(ctx->comm, &rank));
    MPICall(MPI_Comm_size(ctx->comm, &size));
    c->rank = rank;
    c->size = size;

#if defined(DEBUG)
    int min_bs, max_bs;
    MPICall(MPI_Allreduce(&blk, &min_bs, 1, MPI_INT, MPI_MIN, ctx->comm));
    MPICall(MPI_Allreduce(&blk, &max_bs, 1, MPI_INT, MPI_MAX, ctx->comm));
    if (min_bs < 1 || min_bs != max_bs)
        UCFDRaiseError("All ranks must use the same positive block_size");
    if (blk < 0) UCFDRaiseError("Negative block size");
    as_int(product((size_t)bn, (size_t)blk));
#endif

    UCFDInt *range = HOST_ALLOC(UCFDInt, (size_t)size+1);
    build_range(ctx->comm, bn, range);

    bsr->block_elems = product((size_t)blk, (size_t)blk);
    c->blocksize = blk;
    c->vectorized = vectorize;
    c->peer_staging = peer_staging;

    UCFDCall(UCFDSplitBSR(
        bn, blk, rp, ci, range[rank], range[rank+1], bsr
    ));

    /* Host -> Device memory transfer */
    bsr->A.basemat.rowptr = upload_int(bsr->Arp, (size_t)bn + 1);
    bsr->B.basemat.rowptr = upload_int(bsr->Brp, (size_t)bn + 1);
    bsr->d_split_colidx = upload_int(bsr->colidx, (size_t)bsr->nnzb);
    bsr->d_value_source = upload_int(bsr->value_source, (size_t)bsr->nnzb);
    bsr->d_boundary_rows = upload_int(bsr->boundary_rows, (size_t)bsr->n_boundary);
    bsr->d_split_values = device_alloc_real(
        product((size_t)bsr->nnzb, bsr->block_elems)
    );

    bsr->A.bn = bn, bsr->B.bn = bn;
    bsr->A.block = blk, bsr->B.block = blk;
    bsr->A.basemat.colidx = bsr->d_split_colidx;
    bsr->A.basemat.values = bsr->d_split_values;
    bsr->B.basemat.colidx = bsr->d_split_colidx + bsr->nnzb_A;
    bsr->B.basemat.values = bsr->d_split_values + \
                            product((size_t)bsr->nnzb_A, (size_t)bsr->block_elems);

    UCFDCall(UCFDSetCUDAMPIBSRContext(
        &bsr->spmvctx, range, bsr->garray, bsr->n_ghost, blk
    ));
    c->d_send_idx = upload_int(c->send_idx, c->send_total);

    /* Deallocate host arrays */
    free(bsr->Brp);
    free(bsr->value_source);
    free(bsr->boundary_rows);
    free(bsr->garray);
    free(range);

    /* Cached pinned staging supports asynchronous DMA without application-side
     * bounce copies. MPI may read these buffers on the CPU, so avoid write-
     * combined allocations. Actual copy/compute overlap is hardware-dependent. */
    CUDACall(cudaHostAlloc(
        (void **)&c->sbuf, bytes_real(product((size_t)c->send_total, (size_t)blk)+1),
        cudaHostAllocDefault
    ));
    CUDACall(cudaHostAlloc(
        (void **)&c->ghost_x, bytes_real(product((size_t)bsr->n_ghost, (size_t)blk)+1),
        cudaHostAllocDefault
    ));
    c->d_sbuf = device_alloc_real(product((size_t)c->send_total, (size_t)blk));
    c->d_ghost_x = device_alloc_real(product((size_t)bsr->n_ghost, (size_t)blk));
    select_kernels(blk, bsr);

    /* Request setting */
    c->sends = HOST_ALLOC(MPI_Request, c->nsend);
    c->recvs = HOST_ALLOC(MPI_Request, c->nrecv);
    c->completed = HOST_ALLOC(int, c->nrecv);
    for (j=0; j<c->nrecv; ++j)
        MPICall(MPI_Recv_init(
            c->ghost_x + (size_t)c->recv_off[j]*blk,
            scalar_count(c->recv_count[j], blk), MPI_REALTYPE,
            c->recv_nei[j], c->tag, c->comm, &c->recvs[j]
        ));
    for (j = 0; j < c->nsend; ++j)
        MPICall(MPI_Send_init(
            c->sbuf + (size_t)c->send_off[j] * blk,
            scalar_count(c->send_count[j], blk), MPI_REALTYPE,
            c->send_nei[j], c->tag, c->comm, &c->sends[j]
        ));

    /* Three streams: the owned product runs while the copy engines move the
     * halo in both directions. Copy-engine capabilities vary,
     * where supported, inbound and outbound may overlap each other. */
    CUDACall(cudaStreamCreateWithFlags(&c->compute, cudaStreamNonBlocking));
    CUDACall(cudaDeviceGetStreamPriorityRange(&low_priority, &high_priority));
    CUDACall(cudaStreamCreateWithPriority(&c->outbound, cudaStreamNonBlocking, high_priority));
    CUDACall(cudaStreamCreateWithFlags(&c->inbound, cudaStreamNonBlocking));
    CUDACall(cudaEventCreateWithFlags(&c->input_ready, cudaEventDisableTiming));
    CUDACall(cudaEventCreateWithFlags(&c->ghost_ready, cudaEventDisableTiming));
    CUDACall(cudaEventCreateWithFlags(&c->done, cudaEventDisableTiming));
    CUDACall(cudaEventCreateWithFlags(&c->values_ready, cudaEventDisableTiming));
    c->nsend_events = c->peer_staging ? c->nsend : (c->nsend ? 1 : 0);
    c->send_ready = HOST_ALLOC(cudaEvent_t, c->nsend_events);
    for (j = 0; j < c->nsend_events; ++j)
        CUDACall(cudaEventCreateWithFlags(&c->send_ready[j], cudaEventDisableTiming));
    /* ghost_ready is waited on before the first exchange; record it now so the
     * first matvec does not block on an event that was never submitted. */
    CUDACall(cudaEventRecord(c->ghost_ready, c->inbound));
    CUDACall(cudaEventRecord(c->done, c->compute));

    /* `d_va` array must be prepared */
    apply_values(bsr, d_va);

    /* Wait for the initial `d_split_values` array permutation */
    bsr_spmv_synchronize(c);

    m->data             = bsr;
    m->ops->spmv        = SpMV_CUDAMPIBSR;
    m->ops->destroy     = UCFDMPIMatDestroy;

    /* TBU */
    m->ops->update      = UCFDMatUpdateDeviceValues;
    m->ops->cppattern   = UCFDBSRCopyPattern;
    m->ops->cpvalues    = UCFDBSRCopyValues;
    m->ops->cpstream    = UCFDCopyComputeStream;

    UCFDFunctionReturn(UCFD_SUCCESS);
}







