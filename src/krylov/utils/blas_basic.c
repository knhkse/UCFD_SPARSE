#include <omp.h>
#include "blas_basic.h"

void basedaxpy(UCFDInt n, UCFDReal alpha, UCFDReal *x, UCFDReal *y)
{
    UCFDInt i;
    OMPFOR
    for (i=0; i<n; ++i) y[i] += alpha*x[i];
}

UCFDReal basednorm2(MPI_Comm comm, UCFDInt n, UCFDReal *arr)
{
    UCFDInt i;
    UCFDReal sum = 0.0;
    OMPSumReduction(sum)
    for (i=0; i<n; ++i) sum += arr[i]*arr[i];
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_REALTYPE, MPI_SUM, comm);
    return sqrt(sum);
}

UCFDReal baseddot(MPI_Comm comm, UCFDInt n, UCFDReal *x, UCFDReal *y)
{
    UCFDInt i;
    UCFDReal sum = 0.0;
    OMPSumReduction(sum)
    for (i=0; i<n; ++i) sum += x[i]*y[i];
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_REALTYPE, MPI_SUM, comm);
    return sum;
}

void basedscal(UCFDInt n, UCFDReal alpha, UCFDReal *arr)
{
    UCFDInt i;
    OMPFOR
    for (i=0; i<n; ++i) arr[i] = alpha*arr[i];
}

void basedgemvcol(UCFDInt m, UCFDInt n, UCFDInt lda, UCFDReal alpha, UCFDReal *a, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    UCFDInt i, j;
    UCFDReal sum;
    OMPWrapper(j, sum)
    for (i=0; i<m; ++i)
    {
        sum = 0.0;
        for (j=0; j<n; ++j) sum += a[i + j*lda]*x[j];
        y[i] = alpha*sum + beta*y[i];
    }
}

void basedgemvcoltrans(MPI_Comm comm, UCFDInt m, UCFDInt n, UCFDInt lda, UCFDReal *a, UCFDReal *x, UCFDReal *y)
{
    UCFDInt i, j;
    UCFDReal sum;
    OMPWrapper(i, sum)
    for (j=0; j<n; ++j)
    {
        sum = 0.0;
        for (i=0; i<m; ++i) sum += a[i + j*lda]*x[i];
        y[j] = sum;
    }

    MPI_Allreduce(MPI_IN_PLACE, y, (int)n, MPI_REALTYPE, MPI_SUM, comm);
}
