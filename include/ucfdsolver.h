#pragma once

#include "ucfdtypes.h"
#include "ucfdprecon.h"

typedef struct _Solver *Solver;

typedef const char *SolverType;
#define GMRES           "gmres"
#define BICGSTAB        "bicgstab"


#if defined(__cplusplus)
extern "C" {
#endif

UCFD_EXTERN ucfd_status_t UCFDSolverInit(Solver*);
UCFD_EXTERN ucfd_status_t UCFDSolverSetOptions(Solver*, UCFDReal, UCFDReal, UCFDReal, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDSolverDestroy(Solver*);
UCFD_EXTERN ucfd_status_t UCFDSolve(Ctx, Solver, Precon, SpMat, UCFDReal*, UCFDReal*);

UCFD_EXTERN ucfd_status_t UCFDSolverGetResult(Solver, UCFDInt*, UCFDInt*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDSolverTraceResidualHistory(Ctx, Solver);
UCFD_EXTERN ucfd_status_t UCFDSolverGetResidualHistory(Ctx, Solver, UCFDReal*);

UCFD_EXTERN ucfd_status_t UCFDSolverCreateGMRES(Solver*, UCFDInt, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDSolverCreateBICGSTAB(Solver*, UCFDInt);

#if defined(__CUDACC__)
UCFD_EXTERN ucfd_status_t UCFDSolverCreateCUDAGMRES(Solver*, UCFDInt, UCFDInt, UCFDInt, UCFDReal);
#endif

#if defined(__cplusplus)
}
#endif
