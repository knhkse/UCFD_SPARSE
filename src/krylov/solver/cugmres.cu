#include <math.h>
#include "gmres.h"


static ucfd_status_t
arnoldi_cgs2(Solver solver, Precon pc, SpMat A, UCFDInt j, UCFDReal *wnorm)
{
    Solver_CUDAGMRES *ctx = (Solver_CUDAGMRES *)solver->data;
    const UCFDInt k = j + 1;
    const UCFDInt n = ctx->n, m = ctx->restart, ldv = ctx->ldv;
    UCFDReal one=1.0, zero=0.0, neg=-1.0;
    UCFDReal inv;

    UCFDReal *vj        = ctx->d_V + (size_t)(j*ldv);
    UCFDReal *d_h       = ctx->d_proj;            /* offset 0           */
    UCFDReal *d_h2      = ctx->d_proj + (m+1);    /* offset m+1         */
    UCFDReal *d_norm    = ctx->d_proj + 2*(m+1);  /* offset 2(m+1)      */
    UCFDReal *vnext     = ctx->d_V + (size_t)(k*ldv);

    UCFDCall(UCFDSpMV(one, A, vj, zero, ctx->d_w));
    UCFDCall(UCFDPreconApply(pc, ctx->d_w));

    /* ---- CGS2 pass 1:  h = V^T w ;  w -= V h ---- */
    CUBLASCall(cublasDgemv(ctx->handle, CUBLAS_OP_T, n, k, &one,
                           ctx->d_V, ldv, ctx->d_w, 1, &zero, d_h, 1));
    CUBLASCall(cublasDgemv(ctx->handle, CUBLAS_OP_N, n, k, &neg,
                           ctx->d_V, ldv, d_h, 1, &one, ctx->d_w, 1));

    /* ---- CGS2 pass 2 (reorthogonalisation):  h2 = V^T w ;  w -= V h2 ---- */
    CUBLASCall(cublasDgemv(ctx->handle, CUBLAS_OP_T, n, k, &one,
                           ctx->d_V, ldv, ctx->d_w, 1, &zero, d_h2, 1));
    CUBLASCall(cublasDgemv(ctx->handle, CUBLAS_OP_N, n, k, &neg,
                           ctx->d_V, ldv, d_h2, 1, &one, ctx->d_w, 1));

    /* ---- ||w|| written to DEVICE so it rides the single transfer ---- */
    CUBLASCall(cublasSetPointerMode(ctx->handle, CUBLAS_POINTER_MODE_DEVICE));
    CUBLASCall(cublasDnrm2(ctx->handle, n, ctx->d_w, 1, d_norm));
    CUBLASCall(cublasSetPointerMode(ctx->handle, CUBLAS_POINTER_MODE_HOST));

    /* ---- THE single device->host round trip for this iteration ---- */
    CUDACall(cudaMemcpy(ctx->proj_host, ctx->d_proj,
                        (size_t)(2*(m+1)+1)*sizeof(UCFDReal), cudaMemcpyDeviceToHost));

    /* ---- assemble H column j on host:  H[i,j] = h[i] + h2[i] ---- */
    UCFDReal *Hcol = ctx->H + (size_t)j*(m+1);          /* col-major: j*(m+1) */
    for(UCFDInt i=0; i<k; ++i)
        Hcol[i] = ctx->proj_host[i] + ctx->proj_host[(m+1)+i];
    *wnorm = ctx->proj_host[2*(m+1)];
    Hcol[k] = *wnorm;                                  /* H[j+1,j] = ||w||   */

    /* ---- normalise v_{j+1} = w / ||w|| ---- */
    inv = 1.0 / (*wnorm);
    CUBLASCall(cublasDcopy(ctx->handle, n, ctx->d_w, 1, vnext, 1));
    CUBLASCall(cublasDscal(ctx->handle, n, &inv, vnext, 1));

    UCFDFunctionReturn(UCFD_SUCCESS);
}


static void
apply_prev_givens(const UCFDInt m, const UCFDInt j,
                  const UCFDReal *__restrict__ g,
                  UCFDReal *__restrict__ H)
{
    const size_t offset = (size_t)j*(m+1);
    UCFDInt i;
    UCFDReal c, s, h1, h2;
    for (i=0; i<j; ++i)
    {
        c = g[i*2]; s = g[i*2+1];
        h1 = H[i + offset];
        h2 = H[i+1 + offset];
        H[i + offset] = c*h1 + s*h2;
        H[i+1 + offset] = -s*h2 + c*h2;
    }
}

static void
generate_givens(const UCFDInt m, const UCFDInt j,
                UCFDReal *__restrict__ g,
                UCFDReal *__restrict__ H)
{
    const size_t offset = (size_t)j*(m+1);
    UCFDReal h1 = H[j + offset], h2 = H[j+1 + offset];
    UCFDReal rr = hypot(h1, h2);
    UCFDReal c, s;
    if (rr == 0.0) { c = 1.0; s = 0.0; }    /* Degenerate guard */
    else { c = h1/rr; s = h2/rr; }
    g[j*2] = c; g[j*2 + 1] = s;
    H[j + offset] = rr; H[j+1 + offset] = 0.0;
}

static void
update_rhs(const UCFDInt j,
           const UCFDReal *__restrict__ g,
           UCFDReal *__restrict__ y)
{
    const UCFDReal yj = y[j];
    const UCFDReal c = g[j*2], s = g[j*2 + 1];
    y[j] = c*yj;
    y[j+1] = -s*yj;
}

static void
back_substitute(const UCFDInt m, const UCFDInt k,
                const UCFDReal *__restrict__ H,
                UCFDReal *__restrict__ y)
{
    UCFDInt idx, jdx;
    for (idx=k-1; idx>=0; --idx)
    {
        UCFDReal sum = y[idx];
        for (jdx=idx+1; jdx<k; ++jdx)
            sum -= H[idx + (size_t)jdx*(m+1)] * y[jdx];
        y[idx] = sum/H[idx + (size_t)idx*(m+1)];
    }
}

static ucfd_status_t
update_solution(const UCFDInt n, const UCFDInt k,
                Solver_CUDAGMRES *gmres,
                UCFDReal *x)
{
    UCFDReal one = 1.0;

    /* Copy y -> d_y */
    CUDACall(cudaMemcpy(
        gmres->d_y, gmres->y, (size_t)k*sizeof(UCFDReal), cudaMemcpyHostToDevice
    ));

    /* x += Vy */
    CUBLASCall(cublasDgemv(
        gmres->handle, CUBLAS_OP_N, n, k, &one,
        gmres->d_V, ldv, gmres->d_y, 1, &one, x, 1
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t
GMRESSolve(Solver solver, Precon pc, SpMat A, UCFDReal *x, UCFDReal *b)
{
#if defined(DEBUG)
    CheckCUDAPointer(x);
    CheckCUDAPointer(b);
#endif
    UCFDCheckNull(solver->type_name, "Solver must be initialized\n");
    UCFDCheckNull(pc->type_name, "Preconditioner must be initialized\n");
    UCFDCheckNull(A->type_name, "Matrix must be constructed\n");

    Solver_CUDAGMRES *gmres = (Solver_CUDAGMRES *)solver->data;
    CUBLASCall(cublasSetPointerMode(gmres->handle, CUBLAS_POINTER_MODE_HOST));

    const UCFDInt maxiter = solver->maxiter;
    const UCFDReal tol = solver->tol, haptol = solver->haptol;
    const UCFDInt n = gmres->n, m = gmres->restart, ldv = gmres->ldv;
    UCFDReal *H = gmres->H, *g = gmres->g, *y = gmres->y;

    const UCFDInt ld = m + 1;
    UCFDInt iter = 0, j, k;
    UCFDReal wnorm, beta, inv, abeta = 0.0;
    size_t offset;

    while (iter < maxiter)
    {
        /* Compute residual */
        CUBLASCall(cublasDcopy(gmres->handle, n, b, 1, gmres->d_r, 1));
        UCFDCall(UCFDSpMV(-1.0, A, x, 1.0, gmres->d_r));

        /* Convergence check */
        CUBLASCall(cublasDnrm2(gmres->handle, n, gmres->d_r, 1, &abeta));
        solver->ops->record(solver, iter, abeta);
        if (abeta <= tol) {
            solver->stat = CONVERGED;
            break;
        }

        /* Start cycle */
        UCFDCall(UCFDPreconApply(pc, gmres->d_r));
        CUBLASCall(cublasDnrm2(gmres->handle, n, gmres->d_r, 1, &beta));
        gmres->y[0] = beta;
        inv = 1.0/beta;
        CUBLASCall(cublasDcopy(gmres->handle, n, gmres->d_r, 1, gmres->d_V, 1));
        CUBLASCall(cublasDscal(gmres->handle, n, &inv, gmres->d_V, 1));

        k = m;
        for (j=0; j<m; ++j)
        {
            /* Arnoldi iteration */
            UCFDCall(arnoldi_cgs2(solver, pc, A, j, &wnorm));

            /* Givens rotation */
            offset = (size_t)j*ld;

            /* Apply previous givens */
            apply_prev_givens(m, j, g, H);

            /* Generate givens */
            generate_givens(m, j, g, H);

            /* Update rhs */
            update_rhs(j, g, y);

            /* Check convergence */
            if (wnorm < haptol * beta) {
                k = j + 1;
                solver->stat = HAPPYBREAKDOWN;
                break;
            }
        }

        /* Back substitution */
        back_substitute(m, k, H, y);

        /* Update solution */
        UCFDCall(update_solution(n, k, gmres, x));

        iter++;
    }
    if (iter == maxiter) solver->stat = REACH_ITERMAX;
    solver->residual     = abeta;
    solver->itnum        = iter;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t
UCFDDestroyGMRES(Solver solver)
{
    if (!solver) UCFDFunctionReturn(UCFD_SUCCESS);
    Solver_CUDAGMRES *gmres = (Solver_CUDAGMRES *)solver->data;

    free(gmres->H);
    free(gmres->g);
    free(gmres->y);
    free(gmres->proj_host);
    CUDACall(cudaFree(gmres->d_V));
    CUDACall(cudaFree(gmres->d_y));
    CUDACall(cudaFree(gmres->d_w));
    CUDACall(cudaFree(gmres->d_r));
    CUDACall(cudaFree(gmres->d_proj));
    
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline UCFDInt pad_to_16B(UCFDInt n, size_t elem_bytes)
{
    UCFDInt e = (UCFDInt)(16 / elem_bytes);   /* elements per 16 bytes: double->2, float->4 */
    if (e < 1) e = 1;                 /* element already >=16B: no padding needed   */
    return ((n + e - 1) / e) * e;     /* round n up to a multiple of e              */
}

extern "C" ucfd_status_t
UCFDSolverCreateCUDAGMRES(Solver *solver, UCFDInt n, UCFDInt m, UCFDInt maxiter, UCFDReal tol)
{
    UCFDCall(UCFDSolverInit(solver));
    Solver s = *solver;
    s->type_name = GMRES;
    Solver_CUDAGMRES *gmres = (Solver_CUDAGMRES *)calloc(1, sizeof(*gmres));
    UCFDCheckNull(gmres, "GMRES solver allocation failed\n");

#if defined(USE_CUDA)
    CUBLASCall(cublasCreate(&gmres->handle));
#endif

    /* Allocate working arrays */
    /* 1) Host arrays */
    gmres->H            = (UCFDReal*)calloc((size_t)(m+1)*m, sizeof(UCFDReal));
    gmres->g            = (UCFDReal*)calloc((size_t)m*2,    sizeof(UCFDReal));
    gmres->y            = (UCFDReal*)calloc((size_t)(m+1),  sizeof(UCFDReal));
    gmres->proj_host    = (UCFDReal*)calloc((size_t)(2*(m+1)+1),sizeof(UCFDReal));

    /* 2) Device arrays */
    UCFDInt ldv         = pad_to_16B(n, sizeof(UCFDReal));
    CUDACall(cudaMalloc((void**)&gmres->d_V, ldv*(m+1)*sizeof(UCFDReal)));
    CUDACall(cudaMalloc((void**)&gmres->d_y, (m+1)*sizeof(UCFDReal)));
    CUDACall(cudaMalloc((void**)&gmres->d_w, (size_t)n*sizeof(UCFDReal)));
    CUDACall(cudaMalloc((void**)&gmres->d_r, (size_t)n*sizeof(UCFDReal)));
    CUDACall(cudaMalloc((void**)&gmres->d_proj, (size_t)(2*(m+1)+1)*sizeof(UCFDReal)));

    gmres->n                = n;
    gmres->restart          = m;
    gmres->ldv              = ldv;
    s->tol                  = tol;
    s->maxiter              = maxiter;
    s->data                 = gmres;
    s->ops->solve           = GMRESSolve;
    s->ops->destroy         = UCFDDestroyGMRES;

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
