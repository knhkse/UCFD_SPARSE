/**
 * General functions for preconditioner
 */
#include <stdio.h>
#include "ucfdpcimpl.h"
#include "ucfdmatimpl.h"


static void set_diagslots(const UCFDInt n, const UCFDInt *rowptr,
                          const UCFDInt *colidx, UCFDInt *diagslots)
{
    for (UCFDInt i=0; i<n; ++i) {
        UCFDInt st = rowptr[i];
        UCFDInt ed = rowptr[i+1];
        for (UCFDInt j=st; j<ed; ++j) {
            if (i == colidx[j]) {
                diagslots[i] = j;
                continue;
            }
        }
    }
}

ucfd_status_t UCFDPreconCreatefromArrays(Precon *precon, UCFDInt n, UCFDInt nnz, UCFDInt *rowptr, UCFDInt *colidx)
{
    Precon pc = (Precon)calloc(1, sizeof(*pc));
    UCFDCheckNull(pc, "Precon allocation failed\n");

    pc->type_name   = NULL;
    pc->nnz         = nnz;
    pc->rowptr      = rowptr;
    pc->colidx      = colidx;
    pc->diagslots   = malloc((size_t)n * sizeof(UCFDInt));
    set_diagslots(n, pc->rowptr, pc->colidx, pc->diagslots);
    pc->values      = NULL;
    pc->data        = NULL;
    *precon         = pc;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

/**
 * Preconditioner from matrix attributes
 * n : number of (block) rows;
 *      In BSR format, `bn` should be passed
 * nnz : length of `colidx` array
 *      In BSR format, `bnnz` should be passed
 */
ucfd_status_t UCFDPreconCreatefromMatrix(Precon *precon, UCFDInt n, UCFDInt nnz, SpMat mat)
{
    Precon pc = (Precon)calloc(1, sizeof(*pc));
    UCFDCheckNull(pc, "Precon allocation failed\n");

    pc->type_name   = NULL;
    pc->nnz         = nnz;
    mat->ops->cppattern(mat, &pc->rowptr, &pc->colidx);

    /* Construct local-pattern diagslots */
    pc->diagslots       = malloc((size_t)n * sizeof(UCFDInt));
    set_diagslots(n, pc->rowptr, pc->colidx, pc->diagslots);
    pc->values          = NULL;
    pc->data            = NULL;
    *precon             = pc;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPreconPrepare(Precon precon)
{
#if defined(DEBUG)
    UCFDCheckNull(precon->type_name, "Preconditioner type must be set\n");
#endif
    UCFDCall(precon->ops->prepare(precon));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPreconApply(Precon precon, UCFDReal *arr)
{
    UCFDCall(precon->ops->apply(precon, arr));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPreconDestroy(Precon *precon)
{
    if (!precon || !*precon) UCFDFunctionReturn(UCFD_SUCCESS);
    UCFDCall((*precon)->ops->destroy(*precon));
    free((*precon)->data);
    free((*precon)->values);
    free((*precon)->diagslots);
    free(*precon);
    *precon = NULL;
    
    UCFDFunctionReturn(UCFD_SUCCESS);
}
