#include <math.h>
#include "gmres.h"
#include "mpisparsemat.cuh"


static void allreduce_from_device(Solver_CUDAGMRES *gmres, Ctx ctx, const UCFDReal *d_src, UCFDReal *h_dst, int count)
{
    if (gmres->n)
    {
        CUDACall(cudaMemcpyAsync(h_dst, d_src, (size_t)count * sizeof(UCFDReal),
                                cudaMemcpyDeviceToHost, gmres->stream));
        CUDACall(cudaStreamSynchronize(gmres->stream));
    }
    else
        memset(h_dst, 0, (size_t)count * sizeof(UCFDReal));
    if (ctx->size > 1)
        MPICall(MPI_Allreduce(MPI_IN_PLACE, h_dst, count, MPI_DOUBLE, MPI_SUM, ctx->comm));
}

static void push_h2d(Solver_CUDAGMRES *gmres, const UCFDReal *h_src, UCFDReal *d_dst, int count)
{
    if (!gmres->n)
        return;
    CUDACall(cudaMemcpyAsync(d_dst, h_src, (size_t)count * sizeof(UCFDReal),
                            cudaMemcpyHostToDevice, gmres->stream));
}

static void local_sumsq_to_device(Solver_CUDAGMRES *gmres, const UCFDReal *d_v, UCFDReal *d_out)
{
    CUBLASCall(cublasSetPointerMode(gmres->handle, CUBLAS_POINTER_MODE_DEVICE));
    CUBLASCall(cublasDdot(gmres->handle, gmres->n, d_v, 1, d_v, 1, d_out));
    CUBLASCall(cublasSetPointerMode(gmres->handle, CUBLAS_POINTER_MODE_HOST));
}

static double dnrm2_mpi(Solver_CUDAGMRES *gmres, Ctx ctx, const UCFDReal *d_v)
{
    double local = 0.0;
    if (gmres->n)
        CUBLASCall(cublasDdot(gmres->handle, gmres->n, d_v, 1, d_v, 1, &local));
    if (ctx->size > 1)
        MPICall(MPI_Allreduce(MPI_IN_PLACE, &local, 1, MPI_DOUBLE, MPI_SUM, ctx->comm));
    return sqrt(local);
}


static ucfd_status_t
arnoldi_cgs2(Solver solver, Ctx ctx, Precon pc, SpMat A, UCFDInt j, UCFDReal *wnorm)
{
    Solver_CUDAGMRES *gmres = (Solver_CUDAGMRES *)solver->data;
    const UCFDInt n = gmres->n, m = gmres->restart, ldv = gmres->ldv, k = j + 1;
    const int h2_offset = ctx->size == 1 ? k : m + 1;
    UCFDReal one=1.0, zero=0.0, neg=-1.0;

    UCFDReal *vk        = gmres->d_V + (size_t)k * ldv;
    UCFDReal *d_h       = gmres->d_proj;
    UCFDReal *d_h2      = gmres->d_proj + h2_offset;
    UCFDReal *h_h       = gmres->proj_host;
    UCFDReal *h_h2      = gmres->proj_host + h2_offset;
    UCFDReal *Hcol      = gmres->H + (size_t)j * (m + 1);
    UCFDReal hh = 0.0, s, w2sq;
    int i;

    UCFDCall(matrix_spmv(one, A, gmres->d_V + (size_t)j*ldv, zero, vk));
    UCFDCall(apply_precon(pc, vk));

    /* ---- CGS2 pass 1:  h = V^T w ;  w -= V h ---- */
    if (n)
        CUBLASCall(cublasDgemv(gmres->handle, CUBLAS_OP_T, n, k, &one,
                               gmres->d_V, ldv, vk, 1, &zero, d_h, 1));
    if (ctx->size > 1) {
        allreduce_from_device(gmres, ctx, d_h, h_h, k);
        push_h2d(gmres, h_h, d_h, k);
    }

    /* w1 = w - V h : contraction over the REPLICATED dim -> no communication.
     * y is column k, A is columns 0..k-1: disjoint ranges of one allocation. */
    if (n)
        CUBLASCall(cublasDgemv(gmres->handle, CUBLAS_OP_N, n, k, &neg,
                               gmres->d_V, ldv, d_h, 1, &one, vk, 1));

    /* ---- CGS2 pass 2: [h2 | ||w1||^2] = V[:,0..k]^T w1 in ONE GEMV over k+1
            columns — column k IS w1, so its own inner product lands at h2[k].
            The k+1 values are contiguous and ride ONE transfer + ONE reduce --- */
    if (n)
        CUBLASCall(cublasDgemv(gmres->handle, CUBLAS_OP_T, n, k+1, &one,
                               gmres->d_V, ldv, vk, 1, &zero, d_h2, 1));
    if (ctx->size > 1) {
        allreduce_from_device(gmres, ctx, d_h2, h_h2, k+1);
        push_h2d(gmres, h_h2, d_h2, k);
    } else {
        /* Pack [h | h2 | ||w1||^2] without gaps for ONE D2H and ONE host wait.
         * Both device vectors remain valid; no coefficient upload is needed. */
        allreduce_from_device(gmres, ctx, d_h, h_h, 2*k+1);
    }

    /* ---- host, replicated: Hessenberg column + norm without a 3rd reduce --- */
    for (i=0; i<k; ++i)
    {
        Hcol[i] = h_h[i] + h_h2[i];     // CGS2 : column is h + h2
        hh += h_h2[i] * h_h2[i];
    }
    s = h_h2[k];        // Global ||w1||^2
    w2sq = s - hh;

    if (w2sq > 0.25*s)
    {
        /* At most ~2 bits lost: trust it. Because the norm is known BEFORE
         * the correction, v_{j+1} = (w1 - V h2) / ||w2|| is one launch:
         * y = (-1/||w2||) * V h2 + (1/||w2||) * y.  (s > 0 here, so
         * wnorm >= sqrt(DBL_TRUE_MIN) and 1/wnorm cannot overflow.) */
        UCFDReal inv, ninv;
        *wnorm = sqrt(w2sq);
        inv = 1.0 / (*wnorm);
        ninv = -inv;
        if (n)
            CUBLASCall(cublasDgemv(gmres->handle, CUBLAS_OP_N, n, k, &ninv,
                gmres->d_V, ldv, d_h2, 1, &inv, vk, 1));
    } else {
        /* Near-breakdown: the subtraction cancelled, so form w2 unscaled, buy
         * a directly computed norm, and scale afterwards. s and hh are replicated, so
         * every rank takes this branch together. */
        if (n)
            CUBLASCall(cublasDgemv(gmres->handle, CUBLAS_OP_N, n, k, &neg,
                gmres->d_V, ldv, d_h2, 1, &one, vk, 1));
        *wnorm = dnrm2_mpi(gmres, ctx, vk);
        if (n) {
            UCFDReal inv = 1.0 / (*wnorm);
            CUBLASCall(cublasDscal(gmres->handle, n, &inv, vk, 1));
        }
    }
    Hcol[k] = *wnorm;

    UCFDFunctionReturn(UCFD_SUCCESS);
}


static void
apply_prev_givens(const UCFDInt j,
                  const UCFDReal *__restrict__ g,
                  UCFDReal *__restrict__ Hcol)
{
    UCFDInt i;
    UCFDReal c, s, t;
    for (i=0; i<j; ++i)
    {
        c = g[i*2]; s = g[i*2+1];
        t = c*Hcol[i] + s*Hcol[i+1];
        Hcol[i+1] = -s * Hcol[i] + c*Hcol[i+1];
        Hcol[i] = t;
    }
}

static void
givens_generate(const UCFDInt j,
                UCFDReal *__restrict__ g,
                UCFDReal *__restrict__ Hcol)
{
    UCFDReal h1 = Hcol[j], h2 = Hcol[j+1];
    UCFDReal rr = hypot(h1, h2);
    g[j*2] = h1/rr; g[j*2+1] = h2/rr;
    Hcol[j] = g[j*2] * Hcol[j] + g[j*2+1] * Hcol[j+1];
    Hcol[j+1] = 0.0;
}

static void
back_substitute(const UCFDInt k, const UCFDInt ld,
                const UCFDReal *__restrict__ H,
                UCFDReal *__restrict__ y)
{
    for (UCFDInt idx=k-1; idx>=0; --idx) {
        UCFDReal sum = y[idx];
        for (UCFDInt jdx=idx+1; jdx<k; ++jdx)
            sum -= H[idx + jdx*ld] * y[jdx];
        y[idx] = sum/H[idx + idx*ld];
    }
}

static ucfd_status_t
update_solution(const UCFDInt n, const UCFDInt k,
                Solver_CUDAGMRES *gmres,
                UCFDReal *x)
{
    UCFDReal one = 1.0;

    /* Copy y -> d_y */
    CUDACall(cudaMemcpyAsync(
        gmres->d_y, gmres->y, (size_t)k*sizeof(UCFDReal),
        cudaMemcpyHostToDevice, gmres->stream
    ));

    /* x += Vy */
    CUBLASCall(cublasDgemv(
        gmres->handle, CUBLAS_OP_N, n, k, &one,
        gmres->d_V, gmres->ldv, gmres->d_y, 1, &one, x, 1
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t
GMRESSolve(Ctx ctx, Solver solver, Precon pc, SpMat A, UCFDReal *x, UCFDReal *b)
{
#if defined(DEBUG)
    CheckCUDAPointer(x);
    CheckCUDAPointer(b);
    UCFDCheckNull(solver->type_name, "Solver must be initialized\n");
    UCFDCheckNull(pc->type_name, "Preconditioner must be initialized\n");
    UCFDCheckNull(A->type_name, "Matrix must be constructed\n");
#endif
    Solver_CUDAGMRES *gmres = (Solver_CUDAGMRES *)solver->data;
    CUBLASCall(cublasSetPointerMode(gmres->handle, CUBLAS_POINTER_MODE_HOST));

    const UCFDInt n = gmres->n, m = gmres->restart;
    const UCFDInt maxiter = solver->maxiter;
    const UCFDInt maxcycle = (maxiter + m - 1)/m;
    const UCFDInt ld = m + 1;
    
    /* Tolerances */
    const UCFDReal rtol = solver->rtol;
    const UCFDReal atol = solver->atol;
    const UCFDReal dtol = solver->dtol;
    const UCFDReal haptol = solver->haptol;
    
    UCFDReal *H = gmres->H, *g = gmres->g, *y = gmres->y;
    UCFDReal *d_V = gmres->d_V, *d_proj = gmres->d_proj;

    UCFDInt cycle = 0, totit=0, j, k;
    UCFDReal wnorm, inv, beta, c, s, *Hcol;

    UCFDReal beta0 = 0.0;   /* ||M^-1 r_0||, frozen for the whole solve  */
    UCFDReal ttol  = 0.0;   /* the ONE threshold, computed exactly once  */
    UCFDReal rnorm = 0.0;   /* current preconditioned residual estimate  */

    solver->stat = ITERATING;

    /* Initial residual */
    UCFDCall(prepare_precon(pc, A));
    CUBLASCall(cublasDcopy(gmres->handle, n, b, 1, d_V, 1));
    UCFDCall(matrix_spmv(-1.0, A, x, 1.0, d_V));

    while (cycle < maxcycle)
    {
        // beta = start_cycle(gmres, pc, ctx->comm, n, d_r, d_V);
        local_sumsq_to_device(gmres, d_V, d_proj+0);
        if (!strcmp(pc->type_name, NONE)) {
            allreduce_from_device(gmres, ctx, d_proj, gmres->proj_host, 1);
            beta = sqrt(gmres->proj_host[0]);
        } else {
            apply_precon(pc, d_V);
            local_sumsq_to_device(gmres, d_V, d_proj+1);
            allreduce_from_device(gmres, ctx, d_proj, gmres->proj_host, 2);
            beta = sqrt(gmres->proj_host[1]);
        }
        inv = 1.0/beta;
        CUBLASCall(cublasDscal(gmres->handle, n, &inv, d_V, 1));

        y[0] = beta;

        if (cycle == 0) {
            beta0 = beta;
            ttol  = rtol * beta0;
            if (ttol < atol) ttol = atol;

            rnorm = beta;
            solver->ops->record(solver, ctx->rank, 0, rnorm);

            if (rnorm <= ttol) {            /* converged on entry */
                solver->stat = CONVERGED;
                break;
            }
        } else {
            if (fabs(beta - rnorm) > 0.1 * beta0) {
                solver->stat = DIVERGED_BREAKDOWN;
                break;
            }
            rnorm = beta;
            solver->ops->record(solver, ctx->rank, totit, rnorm);
        }

        k = 0;
        for (j=0; j<m && totit < maxiter; ++j)
        {
            /* Arnoldi iteration */
            UCFDCall(arnoldi_cgs2(solver, ctx, pc, A, j, &wnorm));

            /* Givens rotation */
            Hcol = H + (size_t)j * ld;
            apply_prev_givens(j, g, Hcol);
            givens_generate(j, g, Hcol);

            /* Update rhs */
            c = g[j*2]; s = g[j*2+1];
            y[j+1] = -s*y[j];
            y[j] = c*y[j];

            k = j + 1;
            totit++;
            rnorm = fabs(y[j+1]);
            solver->ops->record(solver, ctx->rank, totit, rnorm);

            if (wnorm < haptol) { solver->stat = HAPPYBREAKDOWN; break; }
            if (rnorm < ttol) { solver->stat = CONVERGED; break; }
            if (dtol > 0.0 && rnorm > dtol*beta0) {
                solver->stat = DIVERGED_DTOL;
                break;
            }
        }

        /* Back substitution */
        back_substitute(k, ld, H, y);

        /* Update solution */
        UCFDCall(update_solution(n, k, gmres, x));

        // UCFDCall(calc_residual(solver, A, n, x, b, d_r));
        CUBLASCall(cublasDcopy(gmres->handle, n, b, 1, d_V, 1));
        UCFDCall(matrix_spmv(-1.0, A, x, 1.0, d_V));

        cycle++;
        if (solver->stat != ITERATING) break;
    }
    if (solver->stat == ITERATING) solver->stat = REACH_ITERMAX;
    solver->residual        = rnorm;
    solver->itnum           = totit;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t
UCFDDestroyGMRES(Solver solver)
{
    if (!solver) UCFDFunctionReturn(UCFD_SUCCESS);
    Solver_CUDAGMRES *gmres = (Solver_CUDAGMRES *)solver->data;

    CUDACall(cudaStreamSynchronize(gmres->stream));
    CUDACall(cudaFree(gmres->d_V));
    CUDACall(cudaFree(gmres->d_y));
    CUDACall(cudaFree(gmres->d_proj));
    CUDACall(cudaFreeHost(gmres->proj_host));
    CUDACall(cudaFreeHost(gmres->y));
    free(gmres->H);
    free(gmres->g);

    CUBLASCall(cublasDestroy(gmres->handle));
    
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline UCFDInt pad_to_16B(UCFDInt n, size_t elem_bytes)
{
    UCFDInt e = (UCFDInt)(16 / elem_bytes);   /* elements per 16 bytes: double->2, float->4 */
    if (e < 1) e = 1;                         /* element already >=16B: no padding needed   */
    return ((n + e - 1) / e) * e;             /* round n up to a multiple of e              */
}

extern "C" ucfd_status_t
UCFDSolverCreateCUDAGMRES(Solver *solver, SpMat mat, UCFDInt n, UCFDInt m)
{
    UCFDCall(UCFDSolverInit(solver));
    Solver s = *solver;
    s->type_name = GMRES;
    
    Solver_CUDAGMRES *gmres = (Solver_CUDAGMRES *)calloc(1, sizeof(*gmres));
    UCFDCheckNull(gmres, "GMRES solver allocation failed\n");

    CUBLASCall(cublasCreate(&gmres->handle));
    UCFDCall(UCFDCopyStream(mat, &gmres->stream));
    CUBLASCall(cublasSetStream(gmres->handle, gmres->stream));

    /* Allocate working arrays */
    /* 1) Host arrays */
    gmres->H            = HOST_ALLOC(UCFDReal, (size_t)(m+1)*m);
    gmres->g            = HOST_ALLOC(UCFDReal, (size_t)m*2);
    // Pinned memory for `y` and `proj_host` arrays
    CUDACall(cudaMallocHost((void **)&gmres->y, (size_t)(m+1)*sizeof(UCFDReal)));
    CUDACall(cudaMallocHost((void **)&gmres->proj_host, (size_t)(2*(m+1)+1) * sizeof(UCFDReal)));

    memset(gmres->H, 0, (size_t)(m+1)*m*sizeof(UCFDReal));
    memset(gmres->g, 0, (size_t)m*2*sizeof(UCFDReal));
    memset(gmres->y, 0, (size_t)(m+1)*sizeof(UCFDReal));
    
    /* 2) Device arrays */
    UCFDInt ldv         = pad_to_16B(n, sizeof(UCFDReal));
    CUDACall(cudaMalloc((void**)&gmres->d_V,
        product(product((size_t)ldv, (size_t)(m+1)), sizeof(UCFDReal)) + 1));
    CUDACall(cudaMalloc((void**)&gmres->d_y, (size_t)(m+1)*sizeof(UCFDReal)));
    CUDACall(cudaMalloc((void**)&gmres->d_proj, (size_t)(2*(m+1)+1)*sizeof(UCFDReal)));

    gmres->n                = n;
    gmres->restart          = m;
    gmres->ldv              = ldv;

    s->data                 = gmres;
    s->ops->solve           = GMRESSolve;
    s->ops->destroy         = UCFDDestroyGMRES;
    s->ops->record          = UCFDEmptyKernel;

    // ! Currently, cuBLAS functions are used in default
    s->ops->dcopy           = NULL;
    s->ops->daxpy           = NULL;
    s->ops->dnorm2          = NULL;
    s->ops->ddot            = NULL;
    s->ops->dscal           = NULL;
    s->ops->dgemvcol        = NULL;
    s->ops->dgemvcoltrans   = NULL;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
