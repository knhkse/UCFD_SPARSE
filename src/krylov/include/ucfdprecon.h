#pragma once

#include "ucfd_types.h"
#include "ucfdmat.h"

typedef struct _Precon *Precon;

typedef const char *PreconType;
#define NONE        "none"
#define BILU        "bilu"
#define BLUSGS      "blu-sgs"
#define PBILU       "pbilu"



/* General functions */
// UCFD_EXTERN ucfd_status_t UCFDPreconInitfromMatrix(Precon*, SpMat, UCFDInt*);
UCFD_EXTERN ucfd_status_t UCFDPreconInitfromArrays(Precon*, UCFDInt*, UCFDInt*, UCFDInt*);
UCFD_EXTERN ucfd_status_t UCFDPreconPrep(Precon);
UCFD_EXTERN ucfd_status_t UCFDPreconDestroy(Precon*);

/* Specific functions */
UCFD_EXTERN ucfd_status_t UCFDCreateBILU(Precon*, UCFDInt, UCFDInt, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDCreateBLUSGS(Precon*, UCFDInt, UCFDInt, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDCreatePBILU(Precon*, UCFDInt, UCFDInt, UCFDInt, UCFDInt*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDCreateNonePrecon(Precon*);