#include <math.h>
#include "bicgstab.h"


static ucfd_status_t fused_reduction(UCFDInt n, UCFDReal *restrict rt, UCFDReal *restrict s,
                                     UCFDReal *restrict t, UCFDReal *restrict phi, UCFDReal *restrict rtt,
                                     UCFDReal *restrict ts, UCFDReal *restrict tt, UCFDReal *restrict ss)
{
    UCFDInt k;
    UCFDReal a = 0.0, b = 0.0, c = 0.0, d = 0.0, e = 0.0;
    OMPScheduleStaticSumReduction(a, b, c, d, e)
    for (k=0; k<n; ++k) {
        UCFDReal sk=s[k], tk=t[k], rk=rt[k];
        a += rk * sk;   /* phi = (r~, s) */
        b += rk * tk;   /* (r~, t)       */
        c += tk * sk;   /* (t,  s)       */
        d += tk * tk;   /* (t,  t)       */
        e += sk * sk;   /* (s,  s)       */
    }
    *phi = a; *rtt = b; *ts = c; *tt = d; *ss = e;
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t BICGSTABSolve(Solver solver, Precon pc, SpMat A, UCFDReal *x, UCFDReal *b)
{
    UCFDCheckNull(solver->type_name, "Solver must be initialized\n");
    UCFDCheckNull(pc->type_name, "Preconditioner must be initialized\n");
    UCFDCheckNull(A->type_name, "Matrix must be constructed\n");

    Solver_BICGSTAB *bcs = (Solver_BICGSTAB *)solver->data;
    const UCFDInt maxiter = solver->maxiter, n = bcs->n;
    const UCFDReal tol = solver->tol, haptol = solver->haptol;
    UCFDReal *restrict r = bcs->r, *restrict rt = bcs->rt, *restrict p = bcs->p, \
             *restrict v = bcs->v, *restrict pt = bcs->pt, *restrict s = bcs->s, \
             *restrict shat = bcs->shat, *restrict t = bcs->t;

    UCFDInt iter = 0;
    UCFDReal rho, rhoprev, alpha, beta, omega;
    UCFDReal pi, phi, rtt, ts, tt, ss, rho_new, resnorm=0.0;

    /**
     * Initial residual 
     * 1) r := b
     * 2) r := -A@x + r (b-A@x)
     */
    solver->ops->dcopy(n, r, b);
    UCFDCall(UCFDSpMV(-1.0, A, x, 1.0, r));

    /* rt := r */
    solver->ops->dcopy(n, rt, r);

    /* rho_0 = (rt, r0) */
    rho = solver->ops->ddot(n, rt, r);
    rhoprev = 1.0;
    alpha = 1.0;
    omega = 1.0;

    /* rho breakdown */
    if (fabs(rho) < 1e-30) {
        solver->stat = RHOBREAKDOWN;
        goto done;
    }

    while (iter < maxiter)
    {
        /* direction update p */
        if (iter == 0) solver->ops->dcopy(n, p, r);
        else {
            beta = (rho/rhoprev) * (alpha/omega);
            /* p := r + beta*(p - omega*v) */
            solver->ops->daxpy(n, -omega, v, p);
            solver->ops->dscal(n, beta, p);
            solver->ops->daxpy(n, 1.0, r, p);
        }

        /* pt := inv(M) * p */
        solver->ops->dcopy(n, pt, p);
        UCFDCall(UCFDPreconApply(pc, pt));
        UCFDCall(UCFDSpMV(1.0, A, pt, 0.0, v));

        pi = solver->ops->ddot(n, rt, v);
        if (fabs(pi) < 1e-30) {
            solver->stat = PIBREAKDOWN;
            break;
        }
        alpha = rho/pi;

        /* s := r - alpha*v */
        solver->ops->dcopy(n, s, r);
        solver->ops->daxpy(n, -alpha, v, s);

        /* shat := inv(M)*s */
        solver->ops->dcopy(n, shat, s);
        UCFDCall(UCFDPreconApply(pc, shat));
        UCFDCall(UCFDSpMV(1.0, A, shat, 0.0, t));

        /* Sync 2 : single fused length-5 reduction */
        UCFDCall(fused_reduction(
            n, rt, s, t, &phi, &rtt, &ts, &tt, &ss
        ));
        omega = ts/tt;

        /* solution update : x += alpha*pt + w*s */
        solver->ops->daxpy(n, alpha, pt, x);
        solver->ops->daxpy(n, omega, shat, x);

        /* residual update : r := s - w*t */
        solver->ops->dcopy(n, r, s);
        solver->ops->daxpy(n, -omega, t, r);

        rho_new = phi - omega*rtt;
        resnorm = sqrt(fabs(ss - 2.0*omega*ts + omega*omega*tt));
        solver->ops->record(solver, iter, resnorm);

        if (resnorm < tol) {
            solver->stat = CONVERGED;
            goto done;
        }

        rhoprev = rho;
        rho = rho_new;

        iter++;
    }
done:
    if (iter == maxiter) solver->stat = REACH_ITERMAX;
    solver->itnum = iter;
    solver->residual = resnorm;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t UCFDDestroyBICGSTAB(Solver solver)
{
    if (!solver) UCFDFunctionReturn(UCFD_SUCCESS);
    Solver_BICGSTAB *bcs = (Solver_BICGSTAB *)solver->data;
    free(bcs->r);
    free(bcs->rt);
    free(bcs->p);
    free(bcs->pt);
    free(bcs->v);
    free(bcs->s);
    free(bcs->shat);
    free(bcs->t);

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static struct _SolverOps BICGSTABOps = {
    BICGSTABSolve,
    UCFDDestroyBICGSTAB,
    BLASFUNCS,
    UCFDEmptyKernel
};

ucfd_status_t UCFDSolverCreateBICGSTAB(Solver *solver, UCFDInt n, UCFDInt maxiter, UCFDReal tol)
{
    UCFDCall(UCFDSolverInit(solver));
    Solver s = *solver;
    s->type_name = BICGSTAB;
    Solver_BICGSTAB *bcs = (Solver_BICGSTAB *)calloc(1, sizeof(*bcs));
    UCFDCheckNull(bcs, "BiCGstab solver allocation failed\n");

    bcs->r          = (UCFDReal *)calloc((size_t)n, sizeof(UCFDReal));
    bcs->rt         = (UCFDReal *)calloc((size_t)n, sizeof(UCFDReal));
    bcs->p          = (UCFDReal *)calloc((size_t)n, sizeof(UCFDReal));
    bcs->pt         = (UCFDReal *)calloc((size_t)n, sizeof(UCFDReal));
    bcs->v          = (UCFDReal *)calloc((size_t)n, sizeof(UCFDReal));
    bcs->s          = (UCFDReal *)calloc((size_t)n, sizeof(UCFDReal));
    bcs->shat       = (UCFDReal *)calloc((size_t)n, sizeof(UCFDReal));
    bcs->t          = (UCFDReal *)calloc((size_t)n, sizeof(UCFDReal));
    bcs->n          = n;

    s->tol          = tol;
    s->maxiter      = maxiter;
    s->data         = bcs;
    s->ops[0]       = BICGSTABOps;
    
    UCFDFunctionReturn(UCFD_SUCCESS);
}