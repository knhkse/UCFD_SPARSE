#include "ucfdmatimpl.h"


ucfd_status_t UCFDMatInit(SpMat *mat)
{
    SpMat m = (SpMat)calloc(1, sizeof(*m));
    UCFDNullCheck(m, "Matrix allocation failed\n");

    m->type_name = NULL;
    m->data      = NULL;
    m->n         = 0;
    *mat = m;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDMatMult(UCFDReal alpha, SpMat A, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    UCFDNullCheck(A->type_name, "Matrix type must be set\n");
    UCFDCall(A->ops->spmv(alpha, A, x, beta, y));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDMatDestroy(SpMat *mat)
{
    if (!mat || !*mat) UCFDFunctionReturn(UCFD_SUCCESS);
    UCFDCall((*mat)->ops->destroy(*mat));
    free((*mat)->data);
    free(*mat);
    *mat = NULL;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

