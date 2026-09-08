#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "mpisparsemat.h"
#include "mpihelper.h"


/* Convert a block count to an MPI scalar count without signed overflow. */
static int scalar_count(int nblocks, UCFDInt block_size, MPI_Comm comm)
{
    if (nblocks < 0 || block_size <= 0 || nblocks > INT_MAX / block_size)
    {
        int rank;
        MPI_Comm_rank(comm, &rank);
        fprintf(stderr,
                "rank %d: BSR halo message exceeds MPI int count range\n",
                rank);
        MPI_Abort(comm, MPI_ERR_COUNT);
        abort();
    }
    return nblocks * block_size;
}

static ucfd_status_t UCFDSetBSRMPIContext(UCFDSpMVContext *c,
                                          const UCFDInt *block_range,
                                          const UCFDInt *garray,
                                          UCFDInt ng,
                                          UCFDInt block_size)
{
    int size;
    MPI_Comm_size(c->comm, &size);
    c->nghost = ng;
    c->blocksize = block_size;

    /* All ranks must have same block size */
    int min_block_size;
    int max_block_size;
    MPI_Allreduce(&block_size, &min_block_size, 1, MPI_INT, MPI_MIN, c->comm);
    MPI_Allreduce(&block_size, &max_block_size, 1, MPI_INT, MPI_MAX, c->comm);
    if (min_block_size <= 0 || min_block_size != max_block_size)
    {
        int rank;
        MPI_Comm_rank(c->comm, &rank);
        UCFDRaiseError(
            UCFD_MPI_ERROR, "rank %d: block_size must be positive and equal on all ranks\n"
        );
        MPI_Abort(c->comm, MPI_ERR_ARG);
        abort();
    }

    /* Receive side: unique requested ghost BLOCKS grouped by owning rank. */
    int *request_to = calloc((size_t)size, sizeof(*request_to));
    UCFDInt *request_local =
        malloc((size_t)(ng ? ng : 1) *
               sizeof(*request_local));
    for (UCFDInt k = 0; k < ng; ++k)
    {
        const int p = owner(block_range, size, garray[k]);
        request_local[k] = garray[k] - block_range[p];
        ++request_to[p];
    }

    c->nrecv = 0;
    for (int p = 0; p < size; ++p)
        c->nrecv += (request_to[p] != 0);

    c->recv_nei = malloc((size_t)(c->nrecv ? c->nrecv : 1) *
                              sizeof(*c->recv_nei));
    c->recv_count = malloc((size_t)(c->nrecv ? c->nrecv : 1) *
                              sizeof(*c->recv_count));
    c->recv_off = malloc((size_t)(c->nrecv ? c->nrecv : 1) *
                            sizeof(*c->recv_off));
    for (int p = 0, j = 0, off = 0; p < size; ++p)
    {
        if (request_to[p])
        {
            c->recv_nei[j] = p;
            c->recv_count[j] = request_to[p];
            c->recv_off[j] = off;
            off += request_to[p];
            ++j;
        }
    }

    /* Discover which peers request owned vector BLOCKS from this rank. */
    int *request_from = malloc((size_t)size * sizeof(*request_from));
    MPI_Alltoall(request_to, 1, MPI_INT,
                 request_from, 1, MPI_INT, c->comm);

    c->nsend = 0;
    for (int p = 0; p < size; ++p)
        c->nsend += (request_from[p] != 0);

    c->send_nei = malloc((size_t)(c->nsend ? c->nsend : 1) *
                              sizeof(*c->send_nei));
    c->send_count = malloc((size_t)(c->nsend ? c->nsend : 1) *
                              sizeof(*c->send_count));
    c->send_off = malloc((size_t)(c->nsend ? c->nsend : 1) *
                            sizeof(*c->send_off));
    c->send_total = 0;
    for (int p = 0, j = 0; p < size; ++p)
    {
        if (request_from[p])
        {
            c->send_nei[j] = p;
            c->send_count[j] = request_from[p];
            c->send_off[j] = c->send_total;
            c->send_total += request_from[p];
            ++j;
        }
    }

    c->send_idx = malloc((size_t)(c->send_total ? c->send_total : 1) *
                            sizeof(*c->send_idx));

    /* Exchange owner-local BLOCK indices once during setup. */
    {
        const int nsetup_req = c->nrecv + c->nsend;
        MPI_Request *rq =
            malloc((size_t)(nsetup_req ? nsetup_req : 1) *
                   sizeof(*rq));
        int r = 0;
        for (int j = 0; j < c->nsend; ++j)
            MPI_Irecv(c->send_idx + c->send_off[j],
                      c->send_count[j],
                      MPI_INT, c->send_nei[j], c->tag, c->comm,
                      &rq[r++]);
        for (int j = 0; j < c->nrecv; ++j)
            MPI_Isend(request_local + c->recv_off[j], c->recv_count[j],
                      MPI_INT, c->recv_nei[j], c->tag, c->comm,
                      &rq[r++]);
        if (r) MPI_Waitall(r, rq, MPI_STATUSES_IGNORE);
        free(rq);
    }
    free(request_local);
    free(request_to);
    free(request_from);

    const size_t send_scalars =
        (size_t)c->send_total * (size_t)block_size;
    const size_t ghost_scalars =
        (size_t)ng * (size_t)block_size;
    c->sbuf = malloc((send_scalars ? send_scalars : 1) * sizeof(*c->sbuf));
    c->lvec = malloc((ghost_scalars ? ghost_scalars : 1) *
                          sizeof(*c->lvec));
    c->nreq = c->nrecv + c->nsend;
    c->reqs = malloc((size_t)(c->nreq ? c->nreq : 1) *
                       sizeof(*c->reqs));

    /* Block offsets/counts become scalar offsets/counts at the MPI boundary. */
    for (int j = 0; j < c->nrecv; ++j)
        MPI_Recv_init(c->lvec + (size_t)c->recv_off[j] * (size_t)block_size,
                      scalar_count(c->recv_count[j], block_size, c->comm),
                      MPI_DOUBLE, c->recv_nei[j], c->tag, c->comm,
                      &c->reqs[j]);
    for (int j = 0; j < c->nsend; ++j)
        MPI_Send_init(c->sbuf + (size_t)c->send_off[j] * (size_t)block_size,
                      scalar_count(c->send_count[j], block_size, c->comm),
                      MPI_DOUBLE, c->send_nei[j], c->tag, c->comm,
                      &c->reqs[c->nrecv + j]);

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t UCFDMPIMatDestroy(SpMat mat)
{
    if (!mat) UCFDFunctionReturn(UCFD_SUCCESS);
    MPIBSR *bsr = (MPIBSR *)mat->data;
    UCFDSpMVContext ctx = bsr->spmvctx;

    UCFDCall(UCFDSpMVContextDestroy(&ctx));
    free(bsr->value_dest);
    free(bsr->split_values);
    free(bsr->garray);
    free(bsr->boundary_rows);
    free(bsr->A.basemat.rowptr);
    free(bsr->A.basemat.colidx);
    free(bsr->B.basemat.rowptr);
    free(bsr->B.basemat.colidx);

    UCFDFunctionReturn(UCFD_SUCCESS);
}

/* ---------- phased block-vector halo exchange; no per-call allocation ---------- */
static inline void halo_start(UCFDSpMVContext *c, const UCFDReal *restrict x_local)
{
    const UCFDInt bs = c->blocksize;

    /* Expose receive buffers before local gather work. */
    if (c->nrecv) MPI_Startall(c->nrecv, c->reqs);

    UCFDReal *restrict sbuf = c->sbuf;
    const UCFDInt *restrict send_idx = c->send_idx;
    for (int k = 0; k < c->send_total; ++k)
    {
        const UCFDReal *restrict src =
            x_local + (size_t)send_idx[k] * (size_t)bs;
        UCFDReal *restrict dst = sbuf + (size_t)k * (size_t)bs;
        for (UCFDInt component = 0; component < bs; ++component)
            dst[component] = src[component];
    }

    if (c->nsend) MPI_Startall(c->nsend, c->reqs + c->nrecv);
}


static inline void interior_spmv(UCFDReal alpha,
                                 const BaseBSR *restrict M,
                                 const UCFDReal *restrict x,
                                 UCFDReal beta,
                                 UCFDReal *restrict y)
{
    const size_t bs = (size_t)(M->block);
    const size_t bn = (size_t)(M->bn);
    const size_t block_elems = (size_t)bs * (size_t)bs;
    const UCFDInt *restrict rowptr = M->basemat.rowptr;
    const UCFDInt *restrict colidx = M->basemat.colidx;
    const UCFDReal *restrict val = M->basemat.values;
    UCFDInt ib, row, col, kb;
    UCFDReal accum[bs];

    OMPWrapper(row, col, kb, accum)
    for (ib = 0; ib < bn; ++ib)
    {
        UCFDReal *restrict yblock = y + ib * bs;
        const UCFDInt st = rowptr[ib];
        const UCFDInt end = rowptr[ib + 1];

        for (row = 0; row < bs; ++row)
            accum[row] = 0.0;

        /* Traverse and consume each dense block contiguously. */
        for (kb = st; kb < end; ++kb)
        {
            const UCFDReal *restrict dense =
                val + (size_t)kb * block_elems;
            const UCFDReal *restrict xblock =
                x + (size_t)colidx[kb] * bs;

            for (row = 0; row < bs; ++row)
            {
                const UCFDReal *restrict dense_row =
                    dense + (size_t)row * bs;
                UCFDReal sum = 0.0;

                for (col = 0; col < bs; ++col)
                    sum += dense_row[col] * xblock[col];

                accum[row] += sum;
            }
        }

        for (row = 0; row < bs; ++row)
            yblock[row] = alpha * accum[row] + beta * yblock[row];
    }
}

static inline void boundary_spmv(UCFDReal alpha,
                                 const BaseBSR *restrict M,
                                 const UCFDInt *restrict boundary_rows,
                                 UCFDInt n_boundary,
                                 const UCFDReal *restrict x,
                                 UCFDReal *restrict y)
{
    const size_t bs = (size_t)(M->block);
    const size_t block_elems = (size_t)bs * (size_t)bs;
    const UCFDInt *restrict rowptr = M->basemat.rowptr;
    const UCFDInt *restrict colidx = M->basemat.colidx;
    const UCFDReal *restrict val = M->basemat.values;
    UCFDInt q, kb, row, col;
    UCFDReal accum[bs];

    OMPWrapper(kb, row, col, accum)
    for (q = 0; q < n_boundary; ++q)
    {
        const UCFDInt ib = boundary_rows[q];
        UCFDReal *restrict yblock = y + (size_t)ib * bs;
        const UCFDInt st = rowptr[ib];
        const UCFDInt end = rowptr[ib + 1];

        for (row = 0; row < bs; ++row)
            accum[row] = 0.0;

        /* Traverse and consume each dense ghost block contiguously. */
        for (kb = st; kb < end; ++kb)
        {
            const UCFDReal *restrict dense =
                val + (size_t)kb * block_elems;
            const UCFDReal *restrict xblock =
                x + (size_t)colidx[kb] * bs;

            for (row = 0; row < bs; ++row)
            {
                const UCFDReal *restrict dense_row =
                    dense + (size_t)row * bs;
                UCFDReal sum = 0.0;

                for (col = 0; col < bs; ++col)
                    sum += dense_row[col] * xblock[col];

                accum[row] += sum;
            }
        }

        for (row = 0; row < bs; ++row)
            yblock[row] += alpha * accum[row];
    }
}

static ucfd_status_t SpMV_MPIBSR(UCFDReal alpha, SpMat mat, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    MPIBSR *bsr = (MPIBSR *)mat->data;

    halo_start(&bsr->spmvctx, x);
    interior_spmv(alpha, &bsr->A, x, beta, y);
    halo_wait(&bsr->spmvctx);
    boundary_spmv(alpha, &bsr->B, bsr->boundary_rows, bsr->n_boundary,
                  bsr->spmvctx.lvec, y);
    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t UCFDSplitBSR(UCFDInt bn, UCFDInt blk,
                                  const UCFDInt *rp, const UCFDInt *ci,
                                  const UCFDReal *va, UCFDInt cstart,
                                  UCFDInt cend, MPIBSR *mat)
{
    const UCFDInt nnzb = rp[bn];
    const size_t block_elems = (size_t)blk * (size_t)blk;
    UCFDInt *Arp = malloc((size_t)(bn + 1) * sizeof(*Arp));
    UCFDInt *Brp = malloc((size_t)(bn + 1) * sizeof(*Brp));
    UCFDInt *remote = malloc((size_t)(nnzb ? nnzb : 1) * sizeof(*remote));
    UCFDInt nremote = 0;

    /* Pass 1: count local/remote nonzero BLOCKS and collect ghost blocks. */
    for (UCFDInt ib = 0; ib < bn; ++ib)
    {
        UCFDInt na = 0;
        UCFDInt nb = 0;
        for (UCFDInt kb = rp[ib]; kb < rp[ib + 1]; ++kb)
        {
            const UCFDInt gblock = ci[kb];
            if (gblock >= cstart && gblock < cend)
                ++na;
            else
            {
                ++nb;
                remote[nremote++] = gblock;
            }
        }
        Arp[ib + 1] = na;
        Brp[ib + 1] = nb;
    }

    /* Sorted unique global block columns define compact ghost block slots. */
    qsort(remote, (size_t)nremote, sizeof(*remote), cmp_idx);
    UCFDInt ng = 0;
    for (UCFDInt k = 0; k < nremote; ++k)
        if (k == 0 || remote[k] != remote[k - 1])
            remote[ng++] = remote[k];

    UCFDInt *garray = malloc((size_t)(ng ? ng : 1) * sizeof(*garray));
    if (ng) memcpy(garray, remote, (size_t)ng * sizeof(*garray));
    free(remote);

    UCFDInt n_boundary_block_rows = 0;
    for (UCFDInt ib = 0; ib < bn; ++ib)
        n_boundary_block_rows += (Brp[ib + 1] != 0);

    UCFDInt *boundary_block_rows =
        malloc((size_t)(n_boundary_block_rows ? n_boundary_block_rows : 1) *
               sizeof(*boundary_block_rows));
    for (UCFDInt ib = 0, q = 0; ib < bn; ++ib)
        if (Brp[ib + 1] != 0)
            boundary_block_rows[q++] = ib;

    /* Convert nonzero-block counts to BSR row offsets. */
    Arp[0] = 0;
    Brp[0] = 0;
    for (UCFDInt ib = 0; ib < bn; ++ib)
    {
        Arp[ib + 1] += Arp[ib];
        Brp[ib + 1] += Brp[ib];
    }

    const UCFDInt nnzb_A = Arp[bn];
    const UCFDInt nnzb_B = Brp[bn];
    UCFDInt *Aci = malloc((size_t)(nnzb_A ? nnzb_A : 1) * sizeof(*Aci));
    UCFDInt *Bci = malloc((size_t)(nnzb_B ? nnzb_B : 1) * sizeof(*Bci));
    UCFDInt *value_dest = 
        malloc((size_t)(nnzb ? nnzb : 1)*sizeof(*value_dest));

    /* A.val and B.val are views of one [A blocks | B blocks] allocation. */
    const size_t value_count = (size_t)nnzb * block_elems;
    UCFDReal *split_values =
        malloc((value_count ? value_count : 1) * sizeof(*split_values));
    UCFDReal *Aval = split_values;
    UCFDReal *Bval = split_values + (size_t)nnzb_A * block_elems;

    /*
     * Pass 2: remap block columns and build a permanent original-block to
     * split-block permutation while copying the initial dense blocks.
     */
    UCFDInt ap = 0;
    UCFDInt bp = 0;
    for (UCFDInt ib = 0; ib < bn; ++ib)
    {
        for (UCFDInt kb = rp[ib]; kb < rp[ib + 1]; ++kb)
        {
            const UCFDInt gblock = ci[kb];
            const UCFDReal *src = va + (size_t)kb * block_elems;
            if (gblock >= cstart && gblock < cend)
            {
                Aci[ap] = gblock - cstart;
                value_dest[kb] = ap;
                ++ap;
            }
            else
            {
                Bci[bp] = bsearch_idx(garray, ng, gblock);
                value_dest[kb] = nnzb_A + bp;
                ++bp;
            }
            memcpy(split_values + (size_t)value_dest[kb] * block_elems,
                   src, block_elems * sizeof(*split_values));
        }
    }

    UCFDInt n = bn*blk;
    mat->A = (BaseBSR){{n, Arp, Aci, Aval}, bn, blk};
    mat->B = (BaseBSR){{n, Brp, Bci, Bval}, bn, blk};
    mat->nnzb = nnzb;
    mat->value_dest = value_dest;
    mat->split_values = split_values;
    mat->garray = garray;
    mat->n_ghost = ng;
    mat->boundary_rows = boundary_block_rows;
    mat->n_boundary = n_boundary_block_rows;

    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t UCFDBSRMatUpdate(SpMat mat, UCFDReal *new_values)
{
    MPIBSR *bsr = (MPIBSR *)mat->data;

    const size_t dim2                   = (size_t)bsr->A.block * (size_t)bsr->A.block;
    UCFDReal *restrict dst              = bsr->split_values;
    const UCFDInt *restrict value_dest  = bsr->value_dest;
    const UCFDInt nnzb                  = bsr->nnzb;

    OMPFOR
    for (UCFDInt kb = 0; kb < nnzb; ++kb)
    {
        const UCFDReal *restrict src_block =
            new_values + (size_t)kb * dim2;
        UCFDReal *restrict dst_block =
            dst + (size_t)value_dest[kb] * dim2;

        memcpy(dst_block, src_block, dim2 * sizeof(*dst_block));
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t UCFDBSRCopyPattern(SpMat mat,
                                        UCFDInt **rp_dest,
                                        UCFDInt **ci_dest)
{
    MPIBSR *bsr = (MPIBSR *)mat->data;
    *rp_dest = bsr->A.basemat.rowptr;
    *ci_dest = bsr->A.basemat.colidx;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline ucfd_status_t UCFDBSRCopyValues(SpMat mat, UCFDReal *restrict values)
{
    MPIBSR *bsr = (MPIBSR *)mat->data;
    const UCFDInt val_count = bsr->nnzb * bsr->A.block * bsr->A.block;
    memcpy(values, bsr->A.basemat.values, val_count*sizeof(UCFDReal));

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDMatCreateMPIBSR(SpMat *mat, Ctx ctx,
                                  UCFDInt bn, UCFDInt blk,
                                  UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values)
{
    UCFDCall(UCFDMatInit(mat));
    SpMat m = *mat;
    m->type_name = BSRMPI;

    MPIBSR *bsr = (MPIBSR *)calloc(1, sizeof(*bsr));
    UCFDCheckNull(bsr, "MPIBSR matrix creation failed\n");

    ContextNextTag(ctx, &bsr->spmvctx.tag);
    bsr->spmvctx.comm = ctx->comm;

    /* Prepare : Get range */
    int size=0, rank=0;
    MPI_Comm_size(ctx->comm, &size);
    MPI_Comm_rank(ctx->comm, &rank);
    UCFDInt *range = malloc((size_t)(size + 1)*sizeof(*range));

    bsr->n_local = bn*blk;
    build_range(ctx->comm, bn, range);

    /* Split matrix with interior/boundary region */
    UCFDCall(UCFDSplitBSR(
        bn, blk, rowptr, colidx, values,
        range[rank], range[rank+1], bsr));
    UCFDCall(UCFDSetBSRMPIContext(
        &bsr->spmvctx, range, bsr->garray, bsr->n_ghost, blk
    ));

    free(range);

    m->data             = bsr;
    m->ops->spmv        = SpMV_MPIBSR;
    m->ops->destroy     = UCFDMPIMatDestroy;
    m->ops->update      = UCFDBSRMatUpdate;
    m->ops->cppattern   = UCFDBSRCopyPattern;
    m->ops->cpvalues    = UCFDBSRCopyValues;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
