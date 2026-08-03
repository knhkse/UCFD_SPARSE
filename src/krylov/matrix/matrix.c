#include "ucfdmatimpl.h"


ucfd_status_t UCFDMatInit(SpMat *mat)
{
    SpMat m = (SpMat)calloc(1, sizeof(*m));
    UCFDCheckNull(m, "Matrix allocation failed\n");

    m->type_name    = NULL;
    m->A            = NULL;
    m->B            = NULL;
    *mat            = m;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDMatMult(UCFDReal alpha, SpMat mat, UCFDReal *x, UCFDReal beta, UCFDReal *y)
{
    UCFDCheckNull(mat->type_name, "Matrix type must be set\n");
    UCFDCall(mat->ops->spmv(alpha, mat, x, beta, y));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDMatDestroy(SpMat *mat)
{
    if (!mat || !*mat) UCFDFunctionReturn(UCFD_SUCCESS);
    UCFDCall((*mat)->ops->destroy(*mat));
    free((*mat)->A);
    free((*mat)->B);
    free(*mat);
    *mat = NULL;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

