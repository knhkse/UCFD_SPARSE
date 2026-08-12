#pragma once

#include <mpi.h>
#include "ucfdtypes.h"

typedef struct _Ctx *Ctx;


#if defined(__cplusplus)
extern "C" {
#endif

UCFD_EXTERN ucfd_mpi_t UCFDMPIContextCreate(MPI_Fint, Ctx*);
UCFD_EXTERN ucfd_mpi_t UCFDMPIContextDestroy(Ctx*);

#if defined(__cplusplus)
}
#endif
