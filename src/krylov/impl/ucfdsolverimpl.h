#pragma once

#include "ucfdsolver.h"
#include "ucfdpcimpl.h"
#include "ucfdmatimpl.h"
#include "mpicontext.h"
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
    ucfd_status_t (*solve)(Ctx, Solver, Precon, SpMat, UCFDReal*, UCFDReal*);
    ucfd_status_t (*destroy)(Solver);
    
    /* BLAS functions */
    void (*dcopy)(UCFDInt, UCFDReal*, UCFDReal*);
    void (*daxpy)(UCFDInt, UCFDReal, UCFDReal*, UCFDReal*);
    UCFDReal (*dnorm2)(MPI_Comm, UCFDInt, UCFDReal*);
    UCFDReal (*ddot)(MPI_Comm, UCFDInt, UCFDReal*, UCFDReal*);
    void (*dscal)(UCFDInt, UCFDReal, UCFDReal*);
    void (*dgemvcol)(UCFDInt, UCFDInt, UCFDInt, UCFDReal, UCFDReal*, UCFDReal*, UCFDReal, UCFDReal*);
    void (*dgemvcoltrans)(MPI_Comm, UCFDInt, UCFDInt, UCFDInt, UCFDReal*, UCFDReal*, UCFDReal*);

    /* Trace residual history */
    ucfd_status_t (*record)(Solver, UCFDInt, UCFDInt, UCFDReal);
};

struct _Solver {
    SolverType          type_name;
    UCFDReal            rtol, atol, dtol;
    UCFDReal            haptol;
    UCFDInt             maxiter;
    UCFDInt             itnum;
    UCFDReal            residual, true_residual;
    ucfd_solver_t       stat;
    UCFDReal            *hist_residual;
    void                *data;
    struct _SolverOps   ops[1];
};


/* Encapsulated functions */
static inline ucfd_status_t prepare_precon(Precon precon, SpMat mat)
{
#if defined(DEBUG)
    UCFDCheckNull(precon->type_name, "Preconditioner type must be set\n");
#endif

    // 1) Copy system matrix values -> precon values
    if (precon->values != NULL)
        UCFDCall(mat->ops->cpvalues(mat, precon->values));

    // 2) Prepare preconditioner values
    UCFDCall(precon->ops->prepare(precon));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline ucfd_status_t apply_precon(Precon precon, UCFDReal *__restrict__ b)
{
    UCFDCall(precon->ops->apply(precon, b));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline ucfd_status_t matrix_spmv(UCFDReal alpha, SpMat A, UCFDReal *__restrict__ x,
                                        UCFDReal beta, UCFDReal *__restrict__ y)
{
    UCFDCall(A->ops->spmv(alpha, A, x, beta, y));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

/* Computes r := b - A@x */
static inline ucfd_status_t calc_residual(Solver s, SpMat A, UCFDInt n,
                                          UCFDReal *__restrict__ x, UCFDReal *__restrict__ b, UCFDReal *__restrict__ r)
{
    s->ops->dcopy(n, r, b);
    UCFDCall(matrix_spmv(-1.0, A, x, 1.0, r));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline ucfd_status_t UCFDEmptyKernel(Solver solver, UCFDInt rank, UCFDInt iter, UCFDReal res)
{UCFDFunctionReturn(UCFD_SUCCESS);}
