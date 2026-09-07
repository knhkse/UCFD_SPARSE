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
    UCFDReal            tol;
    UCFDReal            haptol;
    UCFDInt             maxiter;
    UCFDInt             itnum;
    UCFDReal            residual;
    ucfd_solver_t       stat;
    UCFDReal            *hist_residual;
    void                *data;
    struct _SolverOps   ops[1];
};


// ! CUDA -> *restrict (X) *__restrict__ (O) => Need separated function
#if defined(__cplusplus)
extern "C" {
#endif

/* Encapsulated functions */
static inline ucfd_status_t prepare_precon(Precon precon, SpMat mat)
{
#if defined(DEBUG)
    UCFDCheckNull(precon->type_name, "Preconditioner type must be set\n");
#endif

    // 1) Copy system matrix values -> precon values
    UCFDCall(mat->ops->cpvalues(mat, precon->values));

    // 2) Prepare preconditioner
    UCFDCall(precon->ops->prepare(precon));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline ucfd_status_t apply_precon(Precon precon, UCFDReal *restrict b)
{
    UCFDCall(precon->ops->apply(precon, b));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline ucfd_status_t matrix_spmv(UCFDReal alpha, SpMat A, UCFDReal *restrict x,
                                        UCFDReal beta, UCFDReal *restrict y)
{
    UCFDCall(A->ops->spmv(alpha, A, x, beta, y));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

/* Computes r := b - A@x */
static inline ucfd_status_t calc_residual(Solver s, SpMat A, UCFDInt n,
                                          UCFDReal *restrict x, UCFDReal *restrict b, UCFDReal *restrict r)
{
    s->ops->dcopy(n, r, b);
    UCFDCall(matrix_spmv(-1.0, A, x, 1.0, r));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline ucfd_status_t UCFDEmptyKernel(Solver solver, UCFDInt rank, UCFDInt iter, UCFDReal res)
{UCFDFunctionReturn(UCFD_SUCCESS);}

#if defined(__cplusplus)
}
#endif