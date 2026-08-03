#include <math.h>
#include "gmres.h"


static ucfd_status_t
arnoldi_cgs2(Solver solver, Precon pc, SpMat A, UCFDInt j, UCFDReal *wnorm)
{
    Solver_CUDAGMRES *ctx = (Solver_CUDAGMRES *)solver->data;
    const UCFDInt k = j + 1;
    const UCFDInt n = ctx->n;
    const UCFDInt m = ctx->restart, ldv = ctx->ldv;
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
    for(UCFDInt i=0; i<k; i++)
        Hcol[i] = ctx->proj_host[i] + ctx->proj_host[(m+1)+i];
    *wnorm = ctx->proj_host[2*(m+1)];
    Hcol[k] = *wnorm;                                  /* H[j+1,j] = ||w||   */

    /* ---- normalise v_{j+1} = w / ||w|| ---- */
    inv = 1.0 / (*wnorm);
    CUBLASCall(cublasDcopy(ctx->handle, n, ctx->d_w, 1, vnext, 1));
    CUBLASCall(cublasDscal(ctx->handle, n, &inv, vnext, 1));

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t
GMRESSolve(Solver solver, Precon pc, SpMat A, UCFDReal *x, UCFDReal *b)
{
    CheckCUDAPointer(x);
    CheckCUDAPointer(b);

    UCFDCheckNull(solver->type_name, "Solver must be initialized\n");
    UCFDCheckNull(pc->type_name, "Preconditioner must be initialized\n");
    UCFDCheckNull(A->type_name, "Matrix must be constructed\n");
    Solver_CUDAGMRES *gmres = (Solver_CUDAGMRES *)solver->data;
    CUBLASCall(cublasSetPointerMode(gmres->handle, CUBLAS_POINTER_MODE_HOST));

    const UCFDInt n = gmres->n;
    const UCFDInt m = gmres->restart, ldv = gmres->ldv;
    const UCFDInt ld = m + 1;
    UCFDInt iter = 0, i, j, jj, k;
    UCFDReal wnorm, beta, rr, h1, h2, sum, inv, c, s, yj;
    UCFDReal abeta = 0.0, one=1.0;
    size_t offset;

    /* Initial residual */
    CUBLASCall(cublasDcopy(gmres->handle, n, b, 1, gmres->d_r, 1));
    UCFDCall(UCFDSpMV(-1.0, A, x, 1.0, gmres->d_r));

    while (iter < solver->maxiter)
    {
        /* Convergence check */
        CUBLASCall(cublasDnrm2(gmres->handle, n, gmres->d_r, 1, &abeta));
        solver->ops->record(solver, iter, abeta);
        if (abeta <= solver->tol) {
            solver->stat = CONVERGED;
            break;
        }

        /* GMRES step */
        UCFDCall(UCFDPreconApply(pc, gmres->d_r));
        CUBLASCall(cublasDnrm2(gmres->handle, n, gmres->d_r, 1, &beta));
        gmres->y[0] = beta;
        inv = 1.0/beta;
        CUBLASCall(cublasDcopy(gmres->handle, n, gmres->d_r, 1, gmres->d_V, 1));
        CUBLASCall(cublasDscal(gmres->handle, n, &inv, gmres->d_V, 1));

        k = 0;
        for (j=0; j<m; j++)
        {
            /* Arnoldi iteration */
            UCFDCall(arnoldi_cgs2(solver, pc, A, j, &wnorm));

            /* Givens rotation */
            offset = (size_t)j*ld;

            /* Apply previous givens */
            for (i=0; i<j; i++) {
                c = gmres->g[i*2];
                s = gmres->g[i*2+1];
                h1 = gmres->H[i + offset];
                h2 = gmres->H[i+1 + offset];
                gmres->H[i + offset] = c*h1 + s*h2;
                gmres->H[i+1 + offset] = -s*h1 + c*h2;
            }

            /* Generate givens */
            h1 = gmres->H[j + offset];
            h2 = gmres->H[j+1 + offset];
            rr = hypot(h1, h2);
            c = h1/rr; s = h2/rr;
            gmres->g[j*2] = c; gmres->g[j*2+1] = s;
            gmres->H[j + offset] = rr;
            gmres->H[j+1 + offset] = 0.0;

            /* Update rhs */
            yj = gmres->y[j];
            gmres->y[j] = c*yj;
            gmres->y[j+1] = -s*yj;

            k = j + 1;
            if (wnorm < solver->haptol * beta) {
                solver->stat = HAPPYBREAKDOWN;
                break;
            }
        }

        /* Back substitution */
        for (i=k-1; i>=0; --i) {
            sum = gmres->y[i];
            for (jj=i+1; jj<k; ++jj)
                sum -= gmres->H[i + (size_t)jj*ld] * gmres->y[jj];
            gmres->y[i] = sum / gmres->H[i + (size_t)i*ld];
        }

        /* Update solution */
        CUDACall(cudaMemcpy(
            gmres->d_y, gmres->y, (size_t)k*sizeof(UCFDReal), cudaMemcpyHostToDevice
        ));
        CUBLASCall(cublasDgemv(
            gmres->handle, CUBLAS_OP_N, n, k, &one,
            gmres->d_V, ldv, gmres->d_y, 1, &one, x, 1
        ));

        /* Update residual */
        CUBLASCall(cublasDcopy(gmres->handle, n, b, 1, gmres->d_r, 1));
        UCFDCall(UCFDSpMV(-1.0, A, x, 1.0, gmres->d_r));

        iter++;
    }
    if (iter == solver->maxiter) solver->stat = REACH_ITERMAX;
    solver->residual    = abeta;
    solver->itnum     = iter;

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
