#pragma once

#include <string.h>
#include <math.h>
#include "ucfdtypes.h"

static inline void basedcopy(UCFDInt n, UCFDReal *dest, UCFDReal *src)
{
    memcpy(dest, src, sizeof(UCFDReal)*n);
}

void basedaxpy(UCFDInt n, UCFDReal alpha, UCFDReal *x, UCFDReal *y);

UCFDReal basednorm2(MPI_Comm comm, UCFDInt n, UCFDReal *arr);

UCFDReal baseddot(MPI_Comm comm, UCFDInt n, UCFDReal *x, UCFDReal *y);

void basedscal(UCFDInt n, UCFDReal alpha, UCFDReal *arr);

void basedgemvcol(UCFDInt m, UCFDInt n, UCFDInt lda, UCFDReal alpha, UCFDReal *a, UCFDReal *x, UCFDReal beta, UCFDReal *y);

void basedgemvcoltrans(MPI_Comm comm, UCFDInt m, UCFDInt n, UCFDInt lda, UCFDReal *a, UCFDReal *x, UCFDReal *y);







