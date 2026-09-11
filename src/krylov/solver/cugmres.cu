#include <math.h>
#include "gmres.h"


static ucfd_status_t
arnoldi_cgs2(Solver solver, Precon pc, SpMat A, UCFDInt j, UCFDReal *wnorm)
{
    Solver_CUDAGMRES *gmres = (Solver_CUDAGMRES *)solver->data;
    const UCFDInt k = j + 1;
    const UCFDInt n = gmres->n, m = gmres->restart, ldv = gmres->ldv;
    UCFDReal one=1.0, zero=0.0, neg=-1.0;
    UCFDReal inv;

    UCFDReal *vj        = gmres->d_V + (size_t)(j*ldv);
    UCFDReal *d_h       = gmres->d_proj;            /* offset 0           */
    UCFDReal *d_h2      = gmres->d_proj + (m+1);    /* offset m+1         */
    UCFDReal *d_norm    = gmres->d_proj + 2*(m+1);  /* offset 2(m+1)      */
    UCFDReal *vnext     = gmres->d_V + (size_t)(k*ldv);

    UCFDCall(matrix_spmv(one, A, vj, zero, gmres->d_w));
    UCFDCall(apply_precon(pc, gmres->d_w));

    /* ---- CGS2 pass 1:  h = V^T w ;  w -= V h ---- */
    CUBLASCall(cublasDgemv(gmres->handle, CUBLAS_OP_T, n, k, &one,
                           gmres->d_V, ldv, gmres->d_w, 1, &zero, d_h, 1));
    CUBLASCall(cublasDgemv(gmres->handle, CUBLAS_OP_N, n, k, &neg,
                           gmres->d_V, ldv, d_h, 1, &one, gmres->d_w, 1));

    /* ---- CGS2 pass 2 (reorthogonalisation):  h2 = V^T w ;  w -= V h2 ---- */
    CUBLASCall(cublasDgemv(gmres->handle, CUBLAS_OP_T, n, k, &one,
                           gmres->d_V, ldv, gmres->d_w, 1, &zero, d_h2, 1));
    CUBLASCall(cublasDgemv(gmres->handle, CUBLAS_OP_N, n, k, &neg,
                           gmres->d_V, ldv, d_h2, 1, &one, gmres->d_w, 1));

    /* ---- ||w|| written to DEVICE so it rides the single transfer ---- */
    CUBLASCall(cublasSetPointerMode(gmres->handle, CUBLAS_POINTER_MODE_DEVICE));
    CUBLASCall(cublasDnrm2(gmres->handle, n, gmres->d_w, 1, d_norm));
    CUBLASCall(cublasSetPointerMode(gmres->handle, CUBLAS_POINTER_MODE_HOST));

    /* ---- THE single device->host round trip for this iteration ---- */
    CUDACall(cudaMemcpy(gmres->proj_host, gmres->d_proj,
                        (size_t)(2*(m+1)+1)*sizeof(UCFDReal), cudaMemcpyDeviceToHost));

    /* ---- assemble H column j on host:  H[i,j] = h[i] + h2[i] ---- */
    UCFDReal *Hcol = gmres->H + (size_t)j*(m+1);          /* col-major: j*(m+1) */
    for(UCFDInt i=0; i<k; ++i)
        Hcol[i] = gmres->proj_host[i] + gmres->proj_host[(m+1)+i];
    *wnorm = gmres->proj_host[2*(m+1)];
    Hcol[k] = *wnorm;                                  /* H[j+1,j] = ||w||   */

    /* ---- normalise v_{j+1} = w / ||w|| ---- */
    inv = 1.0 / (*wnorm);
    CUBLASCall(cublasDcopy(gmres->handle, n, gmres->d_w, 1, vnext, 1));
    CUBLASCall(cublasDscal(gmres->handle, n, &inv, vnext, 1));

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
    CUDACall(cudaMemcpy(
        gmres->d_y, gmres->y, (size_t)k*sizeof(UCFDReal), cudaMemcpyHostToDevice
    ));

    /* x += Vy */
    CUBLASCall(cublasDgemv(
        gmres->handle, CUBLAS_OP_N, n, k, &one,
        gmres->d_V, gmres->ldv, gmres->d_y, 1, &one, x, 1
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static UCFDReal start_cycle(Solver_CUDAGMRES *gmres, Precon pc, MPI_Comm comm, UCFDInt n,
                            UCFDReal *__restrict__ d_r, UCFDReal *__restrict__ d_V)
{
    UCFDReal beta, inv;

    CUBLASCall(cublasDcopy(gmres->handle, n, d_r, 1, gmres->d_V, 1));
    apply_precon(pc, d_V);
    CUBLASCall(cublasDnrm2(gmres->handle, n, d_V, 1, &beta));
    inv = 1.0/beta;
    CUBLASCall(cublasDscal(gmres->handle, n, &inv, gmres->d_V, 1));

    return beta;
}


extern "C" static ucfd_status_t
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

    const UCFDInt n = gmres->n, m = gmres->restart, ldv = gmres->ldv;
    const UCFDInt maxiter = solver->maxiter;
    const UCFDInt maxcycle = (maxiter + m - 1)/m;
    const UCFDInt ld = m + 1;
    
    /* Tolerances */
    const UCFDReal rtol = solver->rtol;
    const UCFDReal atol = solver->atol;
    const UCFDReal dtol = solver->dtol;
    const UCFDReal haptol = solver->haptol;
    
    UCFDReal *H = gmres->H, *g = gmres->g, *y = gmres->y;
    UCFDReal *d_r = gmres->d_r;

    UCFDInt cycle = 0, totit=0, j, k;
    UCFDReal wnorm, beta, c, s, *Hcol;

    UCFDReal beta0 = 0.0;   /* ||M^-1 r_0||, frozen for the whole solve  */
    UCFDReal ttol  = 0.0;   /* the ONE threshold, computed exactly once  */
    UCFDReal rnorm = 0.0;   /* current preconditioned residual estimate  */

    solver->stat = ITERATING;

    /* Initial residual */
    UCFDCall(prepare_precon(pc, A));
    UCFDCall(calc_residual(solver, A, n, x, b, d_r));

    while (cycle < maxcycle)
    {
        beta = start_cycle(gmres, pc, ctx->comm, n, d_r, gmres->d_V);
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
            UCFDCall(arnoldi_cgs2(solver, pc, A, j, &wnorm));

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

        UCFDCall(calc_residual(solver, A, n, x, b, d_r));

        cycle++;
        if (solver->stat != ITERATING) break;
    }
    if (solver->stat == ITERATING) solver->stat = REACH_ITERMAX;
    solver->residual        = rnorm;
    CUBLASCall(cublasDnrm2(gmres->handle, n, d_r, 1, &solver->true_residual));
    solver->itnum           = totit;

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
UCFDSolverCreateCUDAGMRES(Solver *solver, UCFDInt n, UCFDInt m)
{
    UCFDCall(UCFDSolverInit(solver));
    Solver s = *solver;
    s->type_name = GMRES;

    Solver_CUDAGMRES *gmres = (Solver_CUDAGMRES *)calloc(1, sizeof(*gmres));
    UCFDCheckNull(gmres, "GMRES solver allocation failed\n");

#if defined(__CUDACC__)
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
