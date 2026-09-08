#include "gmres.h"


static ucfd_status_t arnoldi_cgs2(UCFDInt n, UCFDInt m,
                                  Solver solver,
                                  MPI_Comm comm,
                                  UCFDReal *restrict V,
                                  UCFDReal *restrict vj,
                                  UCFDReal *restrict w,
                                  UCFDReal *restrict Hcol,
                                  UCFDReal *restrict h2,
                                  Precon pc, SpMat A, UCFDInt j, UCFDReal *wnorm)
{
    const UCFDInt k = j + 1;
    UCFDReal hsub;

    /* w = inv(M) @ A @ V_j */
    UCFDCall(matrix_spmv(1.0, A, vj, 0.0, w));
    UCFDCall(apply_precon(pc, w));

    /* CGS step 1 */
    solver->ops->dgemvcoltrans(comm, n, k, n, V, w, Hcol);
    solver->ops->dgemvcol(n, k, n, -1.0, V, Hcol, 1.0, w);

    /* CGS step 2 */
    solver->ops->dgemvcoltrans(comm, n, k, n, V, w, h2);
    solver->ops->dgemvcol(n, k, n, -1.0, V, h2, 1.0, w);

    /* H[0:j, j] = h + h2 */
    solver->ops->daxpy(k, 1.0, h2, Hcol);

    /* subdiagonal and new normalized basis vector */
    hsub = solver->ops->dnorm2(comm, n, w);
    Hcol[j + 1] = hsub;
    UCFDReal *vnext = V + (size_t)(j+1) * n;
    solver->ops->dcopy(n, vnext, w);
    solver->ops->dscal(n, 1.0/hsub, vnext);
    *wnorm = hsub;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static void apply_prev_givens(const UCFDInt j,
                              const UCFDReal *restrict cs,
                              const UCFDReal *restrict sn,
                              UCFDReal *restrict Hcol)
{
    UCFDInt i;
    UCFDReal c, s, t;
    for (i=0; i<j; ++i) {
        c = cs[i];
        s = sn[i];
        t = c*Hcol[i] + s*Hcol[i+1];
        Hcol[i+1] = -s * Hcol[i] + c*Hcol[i+1];
        Hcol[i] = t;
    }
}

static void givens_generate(const UCFDInt j,
                            UCFDReal *restrict Hcol,
                            UCFDReal *restrict cs,
                            UCFDReal *restrict sn)
{
    UCFDReal h1 = Hcol[j], h2 = Hcol[j+1];
    UCFDReal rr = hypot(h1, h2);
    cs[j] = h1/rr; sn[j] = h2/rr;
    Hcol[j] = cs[j] * Hcol[j] + sn[j] * Hcol[j+1];
    Hcol[j+1] = 0.0;
}

static void back_substitute(const UCFDInt k, const UCFDInt ld,
                            const UCFDReal *restrict H,
                            UCFDReal *restrict y)
{
    for (UCFDInt i=k-1; i>=0; --i) {
        UCFDReal sum = y[i];
        for (UCFDInt jj=i+1; jj<k; ++jj)
            sum -= H[i+jj*ld] * y[jj];
        y[i] = sum/H[i+i*ld];
    }
}

static UCFDReal start_cycle(Solver solver, Precon pc, MPI_Comm comm, UCFDInt n,
                            UCFDReal *restrict r, UCFDReal *restrict V)
{
    solver->ops->dcopy(n, V, r);
    apply_precon(pc, V);
    UCFDReal beta = solver->ops->dnorm2(comm, n, V);
    solver->ops->dscal(n, 1.0/beta, V);

    return beta;
}


static ucfd_status_t GMRESSolve(Ctx ctx, Solver solver, Precon pc, SpMat A, UCFDReal *x, UCFDReal *b)
{
#if defined(DEBUG)
    UCFDCheckNull(solver->type_name, "Solver must be initialized\n");
    UCFDCheckNull(pc->type_name, "Preconditioner must be initialized\n");
    UCFDCheckNull(A->type_name, "Matrix must be constructed\n");
#endif
    Solver_GMRES *gmres = (Solver_GMRES *)solver->data;

    const UCFDInt n = gmres->n, m = gmres->restart;
    const UCFDInt maxiter = solver->maxiter;
    const UCFDInt maxcycle = gmres->maxcycle;   // maxcycle = (maxiter + m - 1)/m
    const UCFDInt ld = m + 1;

    /* Tolerances */
    const UCFDReal rtol = solver->rtol;
    const UCFDReal atol = solver->atol;
    const UCFDReal dtol = solver->dtol;
    const UCFDReal haptol = solver->haptol;

    UCFDReal *r = gmres->r, *V = gmres->V, *y = gmres->y;
    UCFDReal *cs = gmres->cs, *sn = gmres->sn, *H = gmres->H;

    UCFDInt cycle = 0, totit=0, j, k;
    UCFDReal wnorm, beta, *Hcol;

    UCFDReal beta0 = 0.0;   /* ||M^-1 r_0||, frozen for the whole solve  */
    UCFDReal ttol  = 0.0;   /* the ONE threshold, computed exactly once  */
    UCFDReal rnorm = 0.0;   /* current preconditioned residual estimate  */

    solver->stat = ITERATING;

    /* Initial residual */
    UCFDCall(prepare_precon(pc, A));
    UCFDCall(calc_residual(solver, A, n, x, b, r));

    while (cycle < maxcycle)
    {
        beta = start_cycle(solver, pc, ctx->comm, n, gmres->r, gmres->V);
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
            UCFDCall(arnoldi_cgs2(
                n, m, solver, ctx->comm, gmres->V, gmres->V + j*n, gmres->w, gmres->H + j*ld, gmres->htmp,
                pc, A, j, &wnorm
            ));
            Hcol = H + (size_t)j * ld;

            /* Givens rotation */
            apply_prev_givens(j, cs, sn, Hcol);
            givens_generate(j, Hcol, cs, sn);

            y[j+1] = -sn[j] * y[j];
            y[j] = cs[j] * y[j];

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
        solver->ops->dgemvcol(n, k, n, 1.0, V, y, 1.0, x);

        UCFDCall(calc_residual(solver, A, n, x, b, r));

        cycle++;
        if (solver->stat != ITERATING) break;
    }

    if (solver->stat == ITERATING) solver->stat = REACH_ITERMAX;

    solver->residual        = rnorm;
    solver->true_residual   = solver->ops->dnorm2(ctx->comm, n, r);
    solver->itnum           = totit;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t UCFDDestroyGMRES(Solver solver)
{
    if (!solver) UCFDFunctionReturn(UCFD_SUCCESS);
    Solver_GMRES *gmres = (Solver_GMRES *)solver->data;
    free(gmres->H);
    free(gmres->V);
    free(gmres->y);
    free(gmres->w);
    free(gmres->sn);
    free(gmres->cs);
    free(gmres->htmp);
    free(gmres->r);
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static struct _SolverOps GMRESOps = {
    GMRESSolve,
    UCFDDestroyGMRES,
    BLASFUNCS,
    UCFDEmptyKernel
};

ucfd_status_t UCFDSolverCreateGMRES(Solver *solver, UCFDInt n, UCFDInt m)
{
    UCFDCall(UCFDSolverInit(solver));
    Solver s            = *solver;
    s->type_name        = GMRES;

    Solver_GMRES *gmres = (Solver_GMRES *)calloc(1, sizeof(*gmres));
    UCFDCheckNull(gmres, "GMRES solver allocation failed\n");

    /* Allocate working arrays */
    gmres->H        = (UCFDReal *)calloc((size_t)(m+1)*m, sizeof(UCFDReal));
    gmres->V        = (UCFDReal *)calloc((size_t)n*(m+1), sizeof(UCFDReal));
    gmres->y        = (UCFDReal *)calloc((size_t)(m+1), sizeof(UCFDReal));
    gmres->w        = (UCFDReal *)calloc((size_t)n, sizeof(UCFDReal));
    gmres->sn       = (UCFDReal *)calloc((size_t)(m+1), sizeof(UCFDReal));
    gmres->cs       = (UCFDReal *)calloc((size_t)(m+1), sizeof(UCFDReal));
    gmres->htmp     = (UCFDReal *)calloc((size_t)(m+1), sizeof(UCFDReal));
    gmres->r        = (UCFDReal *)calloc((size_t)n, sizeof(UCFDReal));
    gmres->n        = n;
    gmres->restart  = m;
    gmres->maxcycle = (s->maxiter + m - 1)/m;

    s->data         = gmres;
    s->ops[0]       = GMRESOps;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
