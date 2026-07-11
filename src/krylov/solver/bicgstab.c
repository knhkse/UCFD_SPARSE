#include <math.h>
#include "bicgstab.h"


static ucfd_status_t fused_reduction(UCFDReal *rt, UCFDReal *s, UCFDReal *t,
                                     UCFDInt n, UCFDReal *phi, UCFDReal *rtt,
                                     UCFDReal *ts, UCFDReal *tt, UCFDReal *ss)
{
    UCFDInt k;
    UCFDReal a = 0.0, b = 0.0, c = 0.0, d = 0.0, e = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:a,b,c,d,e)
    for (k=0; k<n; k++) {
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
    UCFDCheckNull(A->type_name, "Matrix must be constructed\n");
    Solver_BICGSTAB *bcs = (Solver_BICGSTAB *)solver->data;
    UCFDInt n = A->n;
    UCFDInt iter = 0;
    UCFDReal rho, rhoprev, alpha, beta, omega;
    UCFDReal pi, phi, rtt, ts, tt, ss, rho_new, resnorm=0.0;

    /**
     * Initial residual 
     * 1) r := b
     * 2) r := -A@x + r (b-A@x)
     */
    solver->ops->dcopy(n, bcs->r, b);
    UCFDCall(UCFDSpMV(-1.0, A, x, 1.0, bcs->r));

    /* rt := r */
    solver->ops->dcopy(n, bcs->rt, bcs->r);

    /* rho_0 = (rt, r0) */
    rho = solver->ops->ddot(n, bcs->rt, bcs->r);
    rhoprev = 1.0;
    alpha = 1.0;
    omega = 1.0;

    /* rho breakdown */
    if (fabs(rho) < 1e-30) {
        solver->stat = RHOBREAKDOWN;
        goto done;
    }

    while (iter < solver->maxiter)
    {
        /* direction update p */
        if (iter == 0) solver->ops->dcopy(n, bcs->p, bcs->r);
        else {
            beta = (rho/rhoprev) * (alpha/omega);
            /* p := r + beta*(p - omega*v) */
            solver->ops->daxpy(n, -omega, bcs->v, bcs->p);
            solver->ops->dscal(n, beta, bcs->p);
            solver->ops->daxpy(n, 1.0, bcs->r, bcs->p);
        }

        /* pt := inv(M) * p */
        solver->ops->dcopy(n, bcs->pt, bcs->p);
        UCFDCall(UCFDPreconApply(pc, bcs->pt));
        UCFDCall(UCFDSpMV(1.0, A, bcs->pt, 0.0, bcs->v));

        pi = solver->ops->ddot(n, bcs->rt, bcs->v);
        if (fabs(pi) < 1e-30) {
            solver->stat = PIBREAKDOWN;
            break;
        }
        alpha = rho/pi;

        /* s := r - alpha*v */
        solver->ops->dcopy(n, bcs->s, bcs->r);
        solver->ops->daxpy(n, -alpha, bcs->v, bcs->s);

        /* shat := inv(M)*s */
        solver->ops->dcopy(n, bcs->shat, bcs->s);
        UCFDCall(UCFDPreconApply(pc, bcs->shat));
        UCFDCall(UCFDSpMV(1.0, A, bcs->shat, 0.0, bcs->t));

        /* Sync 2 : single fused length-5 reduction */
        UCFDCall(fused_reduction(
            bcs->rt, bcs->s, bcs->t, n, &phi, &rtt, &ts, &tt, &ss
        ));
        omega = ts/tt;

        /* solution update : x += alpha*pt + w*s */
        solver->ops->daxpy(n, alpha, bcs->pt, x);
        solver->ops->daxpy(n, omega, bcs->shat, x);

        /* residual update : r := s - w*t */
        solver->ops->dcopy(n, bcs->r, bcs->s);
        solver->ops->daxpy(n, -omega, bcs->t, bcs->r);

        rho_new = phi - omega*rtt;
        resnorm = sqrt(fabs(ss - 2.0*omega*ts + omega*omega*tt));

        if (resnorm < solver->tol) {
            solver->stat = CONVERGED;
            goto done;
        }

        rhoprev = rho;
        rho = rho_new;

        iter++;
    }
done:
    if (iter == solver->maxiter) solver->stat = REACH_ITERMAX;
    solver->subiter = iter;
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
    BLASFUNCS
};

ucfd_status_t UCFDCreateBICGSTAB(Solver *solver, UCFDInt n, UCFDInt maxiter, UCFDReal tol)
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

    s->tol          = tol;
    s->maxiter      = maxiter;
    s->data         = bcs;
    s->ops[0]       = BICGSTABOps;
    
    UCFDFunctionReturn(UCFD_SUCCESS);
}