#pragma once

#include "ucfdtypes.h"
#include "ucfdmpi.h"


typedef struct _SpMat *SpMat;

typedef const char *SpMatType;
#define CSR             "csr"
#define BSR             "bsr"
#define CSRMKL          "mklcsr"
#define BSRMKL          "mklbsr"
#define CSRCUDA         "cudacsr"
#define CSRCUSPARSE     "cusparsecsr"
#define BSRCUDA         "cudabsr"
#define BSRCUSPARSE     "cusparsebsr"
#define CSRMPI          "mpicsr"
#define BSRMPI          "mpibsr"


#if defined(__cplusplus)
extern "C" {
#endif

/* Matrix lifecycle API */
UCFD_EXTERN ucfd_status_t UCFDMatInit(SpMat*);
UCFD_EXTERN ucfd_status_t UCFDMatDestroy(SpMat*);
UCFD_EXTERN ucfd_status_t UCFDMatMult(UCFDReal, SpMat, UCFDReal*, UCFDReal, UCFDReal*);

/* Matrix Setting API */
UCFD_EXTERN ucfd_status_t UCFDMatCreateCSR(SpMat*, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDMatCreateBSR(SpMat*, UCFDInt, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
#if defined(USE_MPI)
UCFD_EXTERN ucfd_status_t UCFDMatCreateMPICSR(Ctx*, SpMat*, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
#endif
#if defined(USE_MKL)
UCFD_EXTERN ucfd_status_t UCFDMatCreateMKLCSR(SpMat*, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDMatCreateMKLBSR(SpMat*, UCFDInt, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
#endif
#if defined(USE_CUDA)
UCFD_EXTERN ucfd_status_t UCFDMatCreateCUDACSR(SpMat*, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDMatCreateCUDABSR(SpMat*, UCFDInt, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
#endif
#if defined(USE_CUSPARSE)
UCFD_EXTERN ucfd_status_t UCFDMatCreateCUSPCSR(SpMat*, UCFDInt, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDMatCreateCUSPBSR(SpMat*, UCFDInt, UCFDInt, UCFDInt, UCFDInt*, UCFDInt*, UCFDReal*);
#endif

#if defined(__cplusplus)
}
#endif
