#pragma once

#include "ucfd_types.h"
#include "ucfdprecon.h"

typedef struct _Solver *Solver;

typedef const char *SolverType;
#define GMRES           "gmres"
#define BICGSTAB        "bicgstab"


#if defined(__cplusplus)
extern "C" {
#endif

UCFD_EXTERN ucfd_status_t UCFDSolverInit(Solver*);
UCFD_EXTERN ucfd_status_t UCFDSolverDestroy(Solver*);
UCFD_EXTERN ucfd_status_t UCFDSolve(Solver, Precon, SpMat, UCFDReal*, UCFDReal*);

UCFD_EXTERN ucfd_status_t UCFDSolverGetResult(Solver, UCFDInt*, UCFDInt*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDSolverTraceResidualHistory(Solver);
UCFD_EXTERN ucfd_status_t UCFDSolverGetResidualHistory(Solver, UCFDReal*);

UCFD_EXTERN ucfd_status_t UCFDSolverCreateGMRES(Solver*, UCFDInt, UCFDInt, UCFDInt, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDSolverCreateBICGSTAB(Solver*, UCFDInt, UCFDInt, UCFDReal);
#if defined(USE_CUDA)
UCFD_EXTERN ucfd_status_t UCFDSolverCreateCUDAGMRES(Solver*, UCFDInt, UCFDInt, UCFDInt, UCFDReal);
#endif

#if defined(__cplusplus)
}
#endif
