#pragma once

#include <mkl.h>
#include "ucfd_types.h"


static inline void mkldcopy(UCFDInt n, UCFDReal *dest, UCFDReal *src)
{
    cblas_dcopy(n, src, 1, dest, 1);
}

static inline void mkldaxpy(UCFDInt n, UCFDReal alpha, UCFDReal *x, UCFDReal *y)
{
    cblas_daxpy(n, alpha, x, 1, y, 1);
}

static inline UCFDReal mkldnorm2(UCFDInt n, UCFDReal *arr)
{
    return cblas_dnrm2(n, arr, 1);
}

static inline UCFDReal mklddot(UCFDInt n, UCFDReal *x, UCFDReal *y)
{
    return cblas_ddot(n, x, 1, y, 1);
}

static inline void mkldscal(UCFDInt n, UCFDReal alpha, UCFDReal *arr)
{
    cblas_dscal(n, alpha, arr, 1);
}

static inline void mkldgemvcol(UCFDInt m, UCFDInt n, UCFDInt lda, UCFDReal alpha, UCFDReal *a, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, alpha, a, lda, x, 1, beta, y, 1);
}

static inline void mkldgemvcoltrans(UCFDInt m, UCFDInt n, UCFDInt lda, UCFDReal alpha, UCFDReal *a, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    cblas_dgemv(CblasColMajor, CblasTrans, m, n, alpha, a, lda, x, 1, beta, y, 1);
}



