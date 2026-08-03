#pragma once

#include "ucfdmat.h"
#include "mpicontext.h"


typedef struct _SpMatOps *SpMatOps;

struct _SpMatOps {
    ucfd_status_t (*spmv)(UCFDReal, SpMat, UCFDReal *, UCFDReal, UCFDReal *);
    ucfd_status_t (*destroy)(SpMat);
};


struct _SpMat {
    SpMatType           type_name;
    void                *A;         // Main matrix
    void                *B;         // Used for MPI communication
    MPIContext          ctx;
    struct _SpMatOps    ops[1];
};

