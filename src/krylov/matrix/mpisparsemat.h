#pragma once
#include "sparsemat.h"
#include "mpicontext.h"


typedef struct {
    BaseCSR         A;
    BaseCSR         B;
    UCFDInt         n_local, n_ghost, n_boundary;
    UCFDInt         *garray;
    UCFDInt         *boundary_rows;
    UCFDMPIContext  ctx;
} MPICSR;




