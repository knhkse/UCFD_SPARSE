#pragma once

#include <mkl.h>
#include "config.h"


static inline void mkldcopy(UCFDInt n, UCFDReal *dest, UCFDReal *src)
{
    cblas_dcopy(n, src, 1, dest, 1);
}

static inline void mkldaxpy(UCFDInt n, UCFDReal alpha, UCFDReal *x, UCFDReal *y)
{
    cblas_daxpy(n, alpha, x, 1, y, 1);
}

static inline UCFDReal mkldnorm2(MPI_Comm comm, UCFDInt n, UCFDReal *arr)
{
    UCFDReal sum = cblas_ddot(n, arr, 1, arr, 1);
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_REALTYPE, MPI_SUM, comm);
    return sqrt(sum);
}

static inline UCFDReal mklddot(MPI_Comm comm, UCFDInt n, UCFDReal *x, UCFDReal *y)
{
    UCFDReal sum = cblas_ddot(n, x, 1, y, 1);
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_REALTYPE, MPI_SUM, comm);
    return sum;
}

static inline void mkldscal(UCFDInt n, UCFDReal alpha, UCFDReal *arr)
{
    cblas_dscal(n, alpha, arr, 1);
}

static inline void mkldgemvcol(UCFDInt m, UCFDInt n, UCFDInt lda, UCFDReal alpha, UCFDReal *a, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, alpha, a, lda, x, 1, beta, y, 1);
}

static inline void mkldgemvcoltrans(MPI_Comm comm, UCFDInt m, UCFDInt n, UCFDInt lda, UCFDReal *a, UCFDReal *x, UCFDReal *y)
{
    cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, a, lda, x, 1, 0.0, y, 1);

    MPI_Allreduce(MPI_IN_PLACE, y, (int)n, MPI_REALTYPE, MPI_SUM, comm);
}



