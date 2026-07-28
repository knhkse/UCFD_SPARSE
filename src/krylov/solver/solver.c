#include <string.h>
#include "ucfdsolverimpl.h"


ucfd_status_t UCFDSolverInit(Solver *solver)
{
    Solver s = (Solver)calloc(1, sizeof(*s));
    UCFDCheckNull(s, "Solver allocation failed\n");

    s->type_name        = NULL;
    s->tol              = 1e-5;
    s->haptol           = 1e-12;
    s->maxiter          = -1;
    s->itnum            = 0;
    s->residual         = 0.0;
    s->stat             = INITIALIZED;
    s->hist_residual    = NULL;
    s->data             = NULL;
    *solver             = s;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDSolverDestroy(Solver *solver)
{
    if (!solver || !*solver) UCFDFunctionReturn(UCFD_SUCCESS);
    UCFDCall((*solver)->ops->destroy(*solver));
    free((*solver)->hist_residual);
    free((*solver)->data);
    free(*solver);
    *solver = NULL;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDSolve(Solver solver, Precon pc, SpMat A, UCFDReal *x, UCFDReal *b)
{
    UCFDCall((solver->ops->solve)(solver, pc, A, x, b));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDSolverGetResult(Solver solver, UCFDInt *stat, UCFDInt *iter, UCFDReal *residual)
{
    *stat       = solver->stat;
    *iter       = solver->itnum;
    *residual   = solver->residual;
    UCFDFunctionReturn(UCFD_SUCCESS);
}

UCFD_INTERN inline ucfd_status_t UCFDSolverStoreResidual(Solver solver, UCFDInt iter, UCFDReal res)
{
    solver->hist_residual[iter] = res;
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDSolverTraceResidualHistory(Solver solver)
{
    solver->hist_residual   = (UCFDReal *)calloc((size_t)solver->maxiter, sizeof(UCFDReal));
    solver->ops->record     = UCFDSolverStoreResidual;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDSolverGetResidualHistory(Solver solver, UCFDReal *hist)
{
    UCFDCheckNull(solver->hist_residual, "Tracing residual unset\n");
    memcpy(hist, solver->hist_residual, solver->maxiter*sizeof(UCFDReal));
    UCFDFunctionReturn(UCFD_SUCCESS);
}
