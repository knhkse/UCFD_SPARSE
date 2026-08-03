#include "gmres.h"


static ucfd_status_t arnoldi_cgs2(UCFDInt n, Solver solver, Precon pc, SpMat A, UCFDInt j, UCFDReal *wnorm)
{
    Solver_GMRES *gmres = (Solver_GMRES *)solver->data;
    const UCFDInt k = j + 1;
    UCFDReal hsub;

    UCFDReal *Vj = gmres->V;
    UCFDReal *vj = gmres->V + (size_t)j * n;
    UCFDReal *w = gmres->w;
    UCFDReal *Hcol = gmres->H + (size_t)j * (gmres->restart + 1);
    UCFDReal *h2 = gmres->htmp;

    /* w = inv(M) @ A @ V_j */
    UCFDCall(UCFDSpMV(1.0, A, vj, 0.0, w));
    UCFDCall(UCFDPreconApply(pc, w));

    /* CGS step 1 */
    solver->ops->dgemvcoltrans(n, k, n, 1.0, Vj, w, 0.0, Hcol);
    solver->ops->dgemvcol(n, k, n, -1.0, Vj, Hcol, 1.0, w);

    /* CGS step 2 */
    solver->ops->dgemvcoltrans(n, k, n, 1.0, Vj, w, 0.0, h2);
    solver->ops->dgemvcol(n, k, n, -1.0, Vj, h2, 1.0, w);

    /* H[0:j, j] = h + h2 */
    solver->ops->daxpy(k, 1.0, h2, Hcol);

    /* subdiagonal and new normalized basis vector */
    hsub = solver->ops->dnorm2(n, w);
    Hcol[j + 1] = hsub;
    UCFDReal *vnext = gmres->V + (size_t)(j+1) * n;
    solver->ops->dcopy(n, vnext, w);
    solver->ops->dscal(n, 1.0/hsub, vnext);
    *wnorm = hsub;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t GMRESSolve(Solver solver, Precon pc, SpMat A, UCFDReal *x, UCFDReal *b)
{
    UCFDCheckNull(solver->type_name, "Solver must be initialized\n");
    UCFDCheckNull(pc->type_name, "Preconditioner must be initialized\n");
    UCFDCheckNull(A->type_name, "Matrix must be constructed\n");

    Solver_GMRES *gmres = (Solver_GMRES *)solver->data;
    const UCFDInt n = gmres->n;
    const UCFDInt m = gmres->restart;
    const UCFDInt ld = m + 1;
    UCFDInt iter = 0, i, j, jj, k;
    UCFDReal wnorm, beta, *Hcol, t, rr, h1, h2, sum, s, c;
    UCFDReal abeta = 0.0; /* Absolute residual */
    
    /**
     * Initial residual 
     * 1) r := b
     * 2) r := -A@x + r (b-A@x)
     */
    solver->ops->dcopy(n, gmres->r, b);
    UCFDCall(UCFDSpMV(-1.0, A, x, 1.0, gmres->r));

    while (iter < solver->maxiter)
    {
        /* Convergence check */
        abeta = solver->ops->dnorm2(n, gmres->r);
        solver->ops->record(solver, iter, abeta);
        if (abeta <= solver->tol) {
            solver->stat = CONVERGED;
            break;
        }

        UCFDCall(UCFDPreconApply(pc, gmres->r));
        beta = solver->ops->dnorm2(n, gmres->r);
        gmres->y[0] = beta;
        solver->ops->dcopy(n, gmres->V, gmres->r);
        solver->ops->dscal(n, 1.0/beta, gmres->V);

        k = 0;
        for (j=0; j<m; j++)
        {
            /* Arnoldi iteration */
            UCFDCall(arnoldi_cgs2(n, solver, pc, A, j, &wnorm));
            Hcol = gmres->H + (size_t)j * ld;

            /* Givens rotation */
            for (i=0; i<j; i++)
            {
                c = gmres->cs[i];
                s = gmres->sn[i];
                t = c * Hcol[i] + s * Hcol[i+1];
                Hcol[i+1] = -s * Hcol[i] + c * Hcol[i+1];
                Hcol[i] = t;
            }
            h1 = Hcol[j]; h2 = Hcol[j+1];
            rr = hypot(h1, h2);
            gmres->cs[j] = h1/rr; gmres->sn[j] = h2/rr;
            Hcol[j] = gmres->cs[j] * Hcol[j] + gmres->sn[j] * Hcol[j+1];
            Hcol[j+1] = 0.0;

            gmres->y[j+1] = -gmres->sn[j] * gmres->y[j];
            gmres->y[j] = gmres->cs[j] * gmres->y[j];

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
                sum -= gmres->H[i+(size_t)jj*ld] * gmres->y[jj];
            gmres->y[i] = sum/gmres->H[i+(size_t)i*ld];
        }
        
        /* Update solution */
        solver->ops->dgemvcol(n, k, n, 1.0, gmres->V, gmres->y, 1.0, x);
        
        /* Update residual */
        solver->ops->dcopy(n, gmres->r, b);
        UCFDCall(UCFDSpMV(-1.0, A, x, 1.0, gmres->r));

        iter++;
    }
    if (iter == solver->maxiter) solver->stat = REACH_ITERMAX;
    solver->residual = abeta;
    solver->itnum = iter;

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

ucfd_status_t UCFDSolverCreateGMRES(Solver *solver, UCFDInt n, UCFDInt m, UCFDInt maxiter, UCFDReal tol)
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

    s->tol          = tol;
    s->maxiter      = maxiter;
    s->data         = gmres;
    s->ops[0]       = GMRESOps;

    UCFDFunctionReturn(UCFD_SUCCESS);
}