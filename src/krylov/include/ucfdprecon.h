#pragma once

#include "ucfdtypes.h"
#include "ucfdmat.h"

typedef struct _Precon *Precon;

typedef const char *PreconType;
#define NONE        "none"
#define ILU         "ilu"
#define BILU        "bilu"
#define BLUSGS      "blu-sgs"
#define PBILU       "pbilu"
#define CUDABILU    "cudabilu"
#define CUDABLUSGS  "cudablusgs"


#if defined(__cplusplus)
extern "C" {
#endif

/* General functions */
UCFD_EXTERN ucfd_status_t UCFDPreconCreatefromArrays(Precon*, UCFDInt*, UCFDInt*, UCFDInt*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDPreconCreateNone(Precon*);
UCFD_EXTERN ucfd_status_t UCFDPreconPrepare(Precon);
UCFD_EXTERN ucfd_status_t UCFDPreconDestroy(Precon*);

/* Specific functions */
UCFD_EXTERN ucfd_status_t UCFDPreconSetILU(Precon *, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDPreconSetBILU(Precon*, UCFDInt, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDPreconSetBLUSGS(Precon*, UCFDInt, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDPreconSetPBILU(Precon*, UCFDInt, UCFDInt, UCFDInt, UCFDInt*);
#if defined(USE_CUDA)
UCFD_EXTERN ucfd_status_t UCFDPreconSetCUDABILU(Precon*, UCFDInt, UCFDInt, UCFDInt, UCFDInt*);
UCFD_EXTERN ucfd_status_t UCFDPreconSetCUDABLUSGS(Precon*, UCFDInt, UCFDInt, UCFDInt, UCFDInt*);
#endif

#if defined(__cplusplus)
}
#endif
