#include <omp.h>
#include "blas_basic.h"

void basedaxpy(UCFDInt n, UCFDReal alpha, UCFDReal *x, UCFDReal *y)
{
    UCFDInt i;
    OMPWrapper()
    for (i=0; i<n; i++) y[i] += alpha*x[i];
}

UCFDReal basednorm2(UCFDInt n, UCFDReal *arr)
{
    UCFDInt i;
    UCFDReal sum = 0.0;
    OMPSumReduction(sum)
    for (i=0; i<n; i++) sum += arr[i]*arr[i];
    return sqrt(sum);
}

UCFDReal baseddot(UCFDInt n, UCFDReal *x, UCFDReal *y)
{
    UCFDInt i;
    UCFDReal sum = 0.0;
    OMPSumReduction(sum)
    for (i=0; i<n; i++) sum += x[i]*y[i];

    return sum;
}

void basedscal(UCFDInt n, UCFDReal alpha, UCFDReal *arr)
{
    UCFDInt i;
    OMPWrapper()
    for (i=0; i<n; i++) arr[i] = alpha*arr[i];
}

void basedgemvcol(UCFDInt m, UCFDInt n, UCFDInt lda, UCFDReal alpha, UCFDReal *a, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    UCFDInt i, j;
    UCFDReal sum;
    OMPWrapper(j, sum)
    for (i=0; i<m; i++)
    {
        sum = 0.0;
        for (j=0; j<n; j++) sum += a[i + j*lda]*x[j];
        y[i] = alpha*sum + beta*y[i];
    }
}

void basedgemvcoltrans(UCFDInt m, UCFDInt n, UCFDInt lda, UCFDReal alpha, UCFDReal *a, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    UCFDInt i, j;
    UCFDReal sum;
    OMPWrapper(i, sum)
    for (j=0; j<n; j++)
    {
        sum = 0.0;
        for (i=0; i<m; i++) sum += a[i + j*lda]*x[i];
        y[j] = alpha*sum + beta*y[j];
    }
}