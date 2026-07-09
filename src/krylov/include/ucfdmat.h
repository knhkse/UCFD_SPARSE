#pragma once

#include "ucfd_types.h"


typedef struct _SpMat *SpMat;

typedef const char *SpMatType;
#define CSR             "csr"
#define BSR             "bsr"
#define MKLCSR          "mklcsr"
#define MKLBSR          "mklbsr"


/* Matrix lifecycle API */
UCFD_EXTERN ucfd_status_t UCFDMatInit(SpMat*);
UCFD_EXTERN ucfd_status_t UCFDMatDestroy(SpMat*);
UCFD_EXTERN ucfd_status_t UCFDMatMult(UCFDReal, SpMat, UCFDReal*, UCFDReal, UCFDReal*);


/* Matrix Setting API */
UCFD_EXTERN ucfd_status_t MatCreateBSR(SpMat*, UCFDInt, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
UCFD_EXTERN ucfd_status_t MatCreateCSR(SpMat*, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
#if defined(USE_MKL)
UCFD_EXTERN ucfd_status_t MatCreateMKLBSR(SpMat*, UCFDInt, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
UCFD_EXTERN ucfd_status_t MatCreateMKLCSR(SpMat*, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
#endif
