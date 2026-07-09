#pragma once

#include <string.h>
#include <math.h>
#include "ucfd_types.h"

static inline void basedcopy(UCFDInt n, UCFDReal *dest, UCFDReal *src)
{
    memcpy(dest, src, sizeof(UCFDReal)*n);
}

void basedaxpy(UCFDInt n, UCFDReal alpha, UCFDReal *x, UCFDReal *y);

UCFDReal basednorm2(UCFDInt n, UCFDReal *arr);

UCFDReal baseddot(UCFDInt n, UCFDReal *x, UCFDReal *y);

void basedscal(UCFDInt n, UCFDReal alpha, UCFDReal *arr);

void basedgemvcol(UCFDInt m, UCFDInt n, UCFDInt lda, UCFDReal alpha, UCFDReal *a, UCFDReal *x, UCFDReal beta, UCFDReal *y);

void basedgemvcoltrans(UCFDInt m, UCFDInt n, UCFDInt lda, UCFDReal alpha, UCFDReal *a, UCFDReal *x, UCFDReal beta, UCFDReal *y);







