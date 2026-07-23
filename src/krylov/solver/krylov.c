#include "ucfdsolverimpl.h"


ucfd_status_t UCFDSolverInit(Solver *solver)
{
    Solver s = (Solver)calloc(1, sizeof(*s));
    UCFDCheckNull(s, "Solver allocation failed\n");

    s->type_name = NULL;
    s->tol       = 1e-5;
    s->haptol    = 1e-12;
    s->maxiter   = 1000;
    s->subiter   = 0;
    s->residual  = 0.0;
    s->stat      = INITIALIZED;
    s->data      = NULL;
    *solver = s;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDSolverDestroy(Solver *solver)
{
    if (!solver || !*solver) UCFDFunctionReturn(UCFD_SUCCESS);
    UCFDCall((*solver)->ops->destroy(*solver));
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

ucfd_status_t UCFDGetSolveResult(Solver solver, UCFDInt *iter, UCFDReal *residual)
{
    *iter = solver->subiter;
    *residual = solver->residual;
    UCFDFunctionReturn(UCFD_SUCCESS);
}