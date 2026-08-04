#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpisparsemat.h"

#define TAG_IDX 1234
#define TAG_VAL 5678


/**
 * Common helper functions
 */
static UCFDInt cmp_idx(const void *a, const void *b)
{
    const UCFDInt x = *(const UCFDInt *)a;
    const UCFDInt y = *(const UCFDInt *)b;
    return (x > y) - (x < y);
}

static UCFDInt bsearch_idx(const UCFDInt *a, UCFDInt m, UCFDInt g)
{
    UCFDInt lo = 0;
    UCFDInt hi = m;

    while (lo < hi)
    {
        const UCFDInt md = lo + (hi - lo) / 2;
        if (a[md] < g)
            lo = md + 1;
        else if (a[md] > g)
            hi = md;
        else
            return md;
    }
    return -1;
}

static size_t owner(const UCFDInt *range, size_t size, UCFDInt g)
{
    size_t lo = 0;
    size_t hi = size;

    while (hi - lo > 1)
    {
        const int md = lo + (hi - lo) / 2;
        if (g < range[md])
            hi = md;
        else
            lo = md;
    }
    return lo;
}

static void build_range(MPI_Comm comm, UCFDInt n_local, UCFDInt *range)
{
    UCFDInt size;

    MPI_Comm_size(comm, &size);
    range[0] = 0;
    MPI_Allgather(&n_local, 1, MPI_INT, range + 1, 1, MPI_INT, comm);
    for (int p = 1; p <= size; ++p)
        range[p] += range[p - 1];
}

/* MPI Context initialization */
static ucfd_status_t UCFDSetMPIContext(MPI_Comm comm,
                                       const UCFDInt *range, const UCFDInt *garray,
                                       UCFDInt n_ghost, UCFDMPIContext *c)
{
    int size, err;

    err = MPI_Comm_dup(comm, &c->comm);
    if (err != MPI_SUCCESS) UCFDFunctionReturn(UCFD_FAILED);
    MPI_Comm_size(c->comm, &size);
    c->nghost = n_ghost;

    /* Receive side: requested counts and owner-local request indices. */
    int *req_to = calloc((size_t)size, sizeof(*req_to));
    UCFDInt *req_local =
        malloc((size_t)(n_ghost ? n_ghost : 1) * sizeof(*req_local));

    for (UCFDInt k = 0; k < n_ghost; ++k)
    {
        const int p = owner(range, size, garray[k]);
        req_local[k] = garray[k] - range[p];
        ++req_to[p];
    }

    c->nrecv = 0;
    for (int p = 0; p < size; ++p) c->nrecv += (req_to[p] != 0);

    c->recv_nei =
        malloc((size_t)(c->nrecv ? c->nrecv : 1) * sizeof(*c->recv_nei));
    c->recv_count =
        malloc((size_t)(c->nrecv ? c->nrecv : 1) * sizeof(*c->recv_count));
    c->recv_off =
        malloc((size_t)(c->nrecv ? c->nrecv : 1) * sizeof(*c->recv_off));

    for (int p = 0, j = 0, off = 0; p < size; ++p)
    {
        if (req_to[p])
        {
            c->recv_nei[j] = p;
            c->recv_count[j] = req_to[p];
            c->recv_off[j] = off;
            off += req_to[p];
            ++j;
        }
    }

    /* Discover send side: who requests from me, and how much */
    int *req_from = malloc((size_t)size * sizeof(*req_from));
    MPI_Alltoall(req_to, 1, MPI_INT, req_from, 1, MPI_INT, c->comm);

    c->nsend = 0;
    for (int p = 0; p < size; ++p)
        c->nsend += (req_from[p] != 0);

    c->send_nei =
        malloc((size_t)(c->nsend ? c->nsend : 1) * sizeof(*c->send_nei));
    c->send_count =
        malloc((size_t)(c->nsend ? c->nsend : 1) * sizeof(*c->send_count));
    c->send_off =
        malloc((size_t)(c->nsend ? c->nsend : 1) * sizeof(*c->send_off));

    c->send_total = 0;
    for (int p = 0, j = 0; p < size; ++p)
    {
        if (req_from[p])
        {
            c->send_nei[j] = p;
            c->send_count[j] = req_from[p];
            c->send_off[j] = c->send_total;
            c->send_total += req_from[p];
            ++j;
        }
    }

    c->send_idx =
        malloc((size_t)(c->send_total ? c->send_total : 1) * sizeof(*c->send_idx));

    /* Exchange index lists once */
    {
        const int nrq = c->nrecv + c->nsend;
        MPI_Request *rq =
            malloc((size_t)(nrq ? nrq : 1) * sizeof(*rq));
        int r = 0;

        for (int j = 0; j < c->nsend; ++j)
            MPI_Irecv(c->send_idx + c->send_off[j],
                      c->send_count[j],
                      MPI_INT,
                      c->send_nei[j],
                      TAG_IDX,
                      c->comm,
                      &rq[r++]);

        for (int j = 0; j < c->nrecv; ++j)
            MPI_Isend(req_local + c->recv_off[j],
                      c->recv_count[j],
                      MPI_INT,
                      c->recv_nei[j],
                      TAG_IDX,
                      c->comm,
                      &rq[r++]);

        if (r)
            MPI_Waitall(r, rq, MPI_STATUSES_IGNORE);
        free(rq);
    }
    free(req_local);
    free(req_to);
    free(req_from);

    c->sbuf =
        malloc((size_t)(c->send_total ? c->send_total : 1) * sizeof(*c->sbuf));
    c->lvec =
        malloc((size_t)(n_ghost ? n_ghost : 1) * sizeof(*c->lvec));

    c->nreq = c->nrecv + c->nsend;
    c->reqs =
        malloc((size_t)(c->nreq ? c->nreq : 1) * sizeof(*c->reqs));

    for (int j = 0; j < c->nrecv; ++j)
        MPI_Recv_init(c->lvec + c->recv_off[j],
                      c->recv_count[j],
                      MPI_DOUBLE,
                      c->recv_nei[j],
                      TAG_VAL,
                      c->comm,
                      &c->reqs[j]);

    for (int j = 0; j < c->nsend; ++j)
        MPI_Send_init(c->sbuf + c->send_off[j],
                      c->send_count[j],
                      MPI_DOUBLE,
                      c->send_nei[j],
                      TAG_VAL,
                      c->comm,
                      &c->reqs[c->nrecv + j]);

    UCFDFunctionReturn(UCFD_SUCCESS);
}


/* Destroy function */
static ucfd_status_t UCFDMPIContextDestroy(UCFDMPIContext *c)
{
    for (int j = 0; j < c->nreq; ++j)
        MPI_Request_free(&c->reqs[j]);
    free(c->reqs);
    free(c->sbuf);
    free(c->lvec);
    free(c->recv_nei);
    free(c->recv_count);
    free(c->recv_off);
    free(c->send_nei);
    free(c->send_count);
    free(c->send_off);
    free(c->send_idx);

    /* Free the private communication context after all requests using it. */
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized && c->comm != MPI_COMM_NULL)
        MPI_Comm_free(&c->comm);

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDMPIMatDestroy(SpMat mat)
{
    if (!mat) UCFDFunctionReturn(UCFD_SUCCESS);
    MPICSR *A = (MPICSR *)mat->data;
    UCFDMPIContext ctx = A->ctx;

    UCFDCall(UCFDMPIContextDestroy(&ctx));
    free(A->garray);
    free(A->boundary_rows);

    UCFDFunctionReturn(UCFD_SUCCESS);
}


/**
 * Halo exchange between processors
 */
static inline void halo_start(UCFDMPIContext *c, const UCFDReal *restrict x_local)
{
    if (c->nrecv) MPI_Startall(c->nrecv, c->reqs);
    UCFDReal *restrict sbuf = c->sbuf;
    const UCFDInt *restrict send_idx = c->send_idx;
    for (int k = 0; k < c->send_total; ++k)
        sbuf[k] = x_local[send_idx[k]];
    if (c->nsend) MPI_Startall(c->nsend, c->reqs + c->nrecv);
}

static inline void halo_wait(UCFDMPIContext *c)
{
    if (c->nreq)
        MPI_Waitall(c->nreq, c->reqs, MPI_STATUSES_IGNORE);
}

/**
 * Local CSR-type SpMV kernel
 */
static inline void interior_spmv(UCFDReal alpha,
                                 const BaseCSR *restrict M,
                                 const UCFDReal *restrict x,
                                 UCFDReal beta,
                                 UCFDReal *restrict y)
{
    const UCFDInt *restrict rowptr = M->rowptr;
    const UCFDInt *restrict colidx = M->colidx;
    const UCFDReal *restrict val = M->values;
    UCFDInt i, k, end;
    UCFDReal sum;

    OMPWrapper(k, end, sum)
    for (i=0; i<M->n; ++i)
    {
        sum = 0.0;
        end = rowptr[i + 1];

        for (k = rowptr[i]; k < end; ++k)
            sum += val[k] * x[colidx[k]];
        y[i] = alpha*sum + beta*y[i];
    }
}

static inline void boundary_spmv(UCFDReal alpha,
                                 const BaseCSR *restrict M,
                                 const UCFDInt *restrict boundary_rows,
                                 UCFDInt n_boundary,
                                 const UCFDReal *restrict x,
                                 UCFDReal *restrict y)
{
    const UCFDInt *restrict rowptr = M->rowptr;
    const UCFDInt *restrict colidx = M->colidx;
    const UCFDReal *restrict val = M->values;
    UCFDInt q, i, k, end;
    UCFDReal sum;

    for (q = 0; q < n_boundary; ++q)
    {
        i = boundary_rows[q];
        sum = 0.0;
        end = rowptr[i + 1];

        for (k = rowptr[i]; k < end; ++k)
            sum += val[k] * x[colidx[k]];
        y[i] += alpha*sum;
    }
}

static ucfd_status_t SpMV_MPICSR(UCFDReal alpha, SpMat mat, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    MPICSR *csr = (MPICSR *)mat->data;

    halo_start(&csr->ctx, x);
    interior_spmv(alpha, &csr->A, x, beta, y);
    halo_wait(&csr->ctx);
    boundary_spmv(alpha, &csr->B, csr->boundary_rows, csr->n_boundary,
                  csr->ctx.lvec, y);
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t UCFDSplitCSR(UCFDInt n, const UCFDInt *rp, const UCFDInt *ci, const UCFDReal *va,
                                  UCFDInt cstart, UCFDInt cend, MPICSR *mat)
{
    const UCFDInt nnz = rp[n];
    UCFDInt *Arp = malloc((size_t)(n + 1) * sizeof(*Arp));
    UCFDInt *Brp = malloc((size_t)(n + 1) * sizeof(*Brp));
    UCFDInt *rem = malloc((size_t)(nnz ? nnz : 1) * sizeof(*rem));
    UCFDInt nrem = 0;

    /* Pass 1: count per-row A/B nonzeros and collect remote columns */
    for (UCFDInt i = 0; i < n; ++i)
    {
        UCFDInt na = 0;
        UCFDInt nb = 0;

        for (UCFDInt k = rp[i]; k < rp[i + 1]; ++k)
        {
            const UCFDInt g = ci[k];
            if (g >= cstart && g < cend)
                ++na;
            else
            {
                ++nb;
                rem[nrem++] = g;
            }
        }
        Arp[i + 1] = na;
        Brp[i + 1] = nb;
    }

    /* Build garray by sorting and uniquing the remote global columns. */
    qsort(rem, (size_t)nrem, sizeof(*rem), cmp_idx);
    UCFDInt ng = 0;
    for (UCFDInt k = 0; k < nrem; ++k)
        if (k == 0 || rem[k] != rem[k - 1])
            rem[ng++] = rem[k];

    UCFDInt *garray = malloc((size_t)(ng ? ng : 1) * sizeof(*garray));
    if (ng) memcpy(garray, rem, (size_t)ng * sizeof(*garray));
    free(rem);

    /* Build boundary rows */
    UCFDInt n_boundary_rows = 0;
    for (UCFDInt i = 0; i < n; ++i)
        n_boundary_rows += (Brp[i + 1] != 0);

    UCFDInt *boundary_rows =
        malloc((size_t)(n_boundary_rows ? n_boundary_rows : 1) * sizeof(*boundary_rows));
    for (UCFDInt i = 0, j = 0; i < n; ++i)
        if (Brp[i + 1] != 0)
            boundary_rows[j++] = i;

    /* Convert per-row counts to CSR prefix sums. */
    Arp[0] = 0;
    Brp[0] = 0;
    for (UCFDInt i = 0; i < n; ++i)
    {
        Arp[i + 1] += Arp[i];
        Brp[i + 1] += Brp[i];
    }

    const UCFDInt nA = Arp[n];
    const UCFDInt nB = Brp[n];
    UCFDInt *Aci = malloc((size_t)(nA ? nA : 1) * sizeof(*Aci));
    UCFDInt *Bci = malloc((size_t)(nB ? nB : 1) * sizeof(*Bci));
    UCFDReal *Av = malloc((size_t)(nA ? nA : 1) * sizeof(*Av));
    UCFDReal *Bv = malloc((size_t)(nB ? nB : 1) * sizeof(*Bv));

    /* Pass 2: fill A and B, remapping global columns to local slots */
    UCFDInt ap = 0;
    UCFDInt bp = 0;
    for (UCFDInt i = 0; i < n; ++i)
    {
        for (UCFDInt k = rp[i]; k < rp[i + 1]; ++k)
        {
            const UCFDInt g = ci[k];
            if (g >= cstart && g < cend)
            {
                Aci[ap] = g - cstart;
                Av[ap] = va[k];
                ++ap;
            }
            else
            {
                Bci[bp] = bsearch_idx(garray, ng, g);
                Bv[bp] = va[k];
                ++bp;
            }
        }
    }

    mat->A = (BaseCSR){n, Arp, Aci, Av};
    mat->B = (BaseCSR){n, Brp, Bci, Bv};
    mat->garray = garray;
    mat->n_ghost = ng;
    mat->boundary_rows = boundary_rows;
    mat->n_boundary = n_boundary_rows;

    UCFDFunctionReturn(UCFD_SUCCESS);
}


/**
 * CSR MPI Matrix
 */
ucfd_status_t UCFDMatCreateMPICSR(MPI_Fint fcomm, SpMat *mat, UCFDInt n, UCFDInt *rowptr, UCFDInt *colidx, UCFDReal *values)
{
    /* Check MPI preparation */
    int flag = 0;
    MPI_Comm comm_in;

    // TODO : Apply `mpi_spmv_ctypes.c` mechanism
    MPI_Initialized(&flag);
    if (!flag) UCFDFunctionReturn(UCFD_FAILED);
    comm_in = MPI_Comm_f2c(fcomm);
    if (comm_in == MPI_COMM_NULL) UCFDFunctionReturn(UCFD_FAILED);

    /* Initialize matrix object */
    UCFDCall(UCFDMatInit(mat));
    SpMat m = *mat;
    m->type_name = CSRMPI;

    MPICSR *csr = (MPICSR *)calloc(1, sizeof(*csr));
    UCFDCheckNull(csr, "MPICSR matrix creation failed\n");

    /* Prepare : Get range */
    int size=0, rank=0;
    MPI_Comm_size(comm_in, &size);
    MPI_Comm_rank(comm_in, &rank);
    UCFDInt *range = malloc((size_t)(size + 1)*sizeof(*range));

    csr->n_local = n;
    build_range(comm_in, n, range);

    /* Split matrix with interior/boundary region */
    UCFDCall(UCFDSplitCSR(
        n, rowptr, colidx, values,
        range[rank], range[rank+1], csr));
    UCFDCall(UCFDSetMPIContext(
        comm_in, range,
        csr->garray, csr->n_ghost, &csr->ctx
    ));

    free(range);

    m->data         = csr;
    m->ops->spmv    = SpMV_MPICSR;
    m->ops->destroy = UCFDMPIMatDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
