#pragma once

#include "ucfdmat.h"


typedef struct _SpMatOps *SpMatOps;
struct _SpMatOps {
    ucfd_status_t (*spmv)(UCFDReal, SpMat, UCFDReal *, UCFDReal, UCFDReal *);
    ucfd_status_t (*destroy)(SpMat);
};


struct _SpMat {
    SpMatType   type_name;
    UCFDInt     initialized;

    /* Total matrix size to allocate working arrays */
    UCFDInt     n;

    /* Sparse arrays for CSR-based format */
    UCFDInt     *rowptr;
    UCFDInt     *colidx;
    UCFDReal    *values;
    void        *data;
    struct _SpMatOps ops[1];
};






