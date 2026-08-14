#pragma once

#include "ucfdsolver.h"
#include "ucfdpcimpl.h"
#include "ucfdmatimpl.h"
#include "blas_basic.h"

#if defined(USE_MKL)
    #include "blas_mkl.h"
    #define BLASFUNCS       \
        mkldcopy,           \
        mkldaxpy,           \
        mkldnorm2,          \
        mklddot,            \
        mkldscal,           \
        mkldgemvcol,        \
        mkldgemvcoltrans
#else
    #define BLASFUNCS       \
        basedcopy,          \
        basedaxpy,          \
        basednorm2,         \
        baseddot,           \
        basedscal,          \
        basedgemvcol,       \
        basedgemvcoltrans
#endif

typedef struct _SolverOps *SolverOps;
struct _SolverOps {
    /* Solver functions */
    ucfd_status_t (*solve)(Solver, Precon, SpMat, UCFDReal*, UCFDReal*);
    ucfd_status_t (*destroy)(Solver);
    
    /* BLAS functions */
    void (*dcopy)(UCFDInt, UCFDReal*, UCFDReal*);
    void (*daxpy)(UCFDInt, UCFDReal, UCFDReal*, UCFDReal*);
    UCFDReal (*dnorm2)(UCFDInt, UCFDReal*);
    UCFDReal (*ddot)(UCFDInt, UCFDReal*, UCFDReal*);
    void (*dscal)(UCFDInt, UCFDReal, UCFDReal*);
    void (*dgemvcol)(UCFDInt, UCFDInt, UCFDInt, UCFDReal, UCFDReal*, UCFDReal*, UCFDReal, UCFDReal*);
    void (*dgemvcoltrans)(UCFDInt, UCFDInt, UCFDInt, UCFDReal, UCFDReal*, UCFDReal*, UCFDReal, UCFDReal*);

    /* et cetera... */
    ucfd_status_t (*record)(Solver, UCFDInt, UCFDReal);
};

struct _Solver {
    SolverType type_name;
    UCFDReal tol;
    UCFDReal haptol;
    UCFDInt maxiter;
    UCFDInt itnum;
    UCFDReal residual;
    ucfd_solver_t stat;
    UCFDReal *hist_residual;
    void *data;
    struct _SolverOps ops[1];
};


/* Encapsulated functions */
static inline ucfd_status_t UCFDPreconApply(Precon precon, UCFDReal *b)
{
    UCFDCall(precon->ops->apply(precon, b));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline ucfd_status_t UCFDSpMV(UCFDReal alpha, SpMat A, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    UCFDCall(A->ops->spmv(alpha, A, x, beta, y));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

/* Computes r := b - A@x */
static inline ucfd_status_t UCFDCalcResidual(Solver s, SpMat A, UCFDInt n, UCFDReal *x, UCFDReal *b, UCFDReal *r)
{
    s->ops->dcopy(n, r, b);
    UCFDCall(UCFDSpMV(-1.0, A, x, 1.0, r));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline ucfd_status_t UCFDEmptyKernel(Solver solver, UCFDInt iter, UCFDReal res)
{UCFDFunctionReturn(UCFD_SUCCESS);}
