#pragma once

#include "ucfdmat.h"


typedef struct _SpMatOps *SpMatOps;

struct _SpMatOps {
    ucfd_status_t (*spmv)(UCFDReal, SpMat, UCFDReal *, UCFDReal, UCFDReal *);
    ucfd_status_t (*destroy)(SpMat);
    ucfd_status_t (*update)(SpMat, UCFDReal*);
    ucfd_status_t (*cppattern)(SpMat, UCFDInt**, UCFDInt**, UCFDReal**);
    ucfd_status_t (*cpvalues)(SpMat, UCFDReal*);
};


struct _SpMat {
    SpMatType           type_name;
    void                *data;
    struct _SpMatOps    ops[1];
};


#if defined(__cplusplus)
extern "C" {
#endif

UCFD_INTERN ucfd_status_t UCFDMatCopyPattern(SpMat, UCFDInt**, UCFDInt**, UCFDReal**);

#if defined(__cplusplus)
}
#endif
