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

ucfd_status_t UCFDPreconCreatefromArrays(Precon *precon, UCFDInt *rowptr, UCFDInt *colidx, UCFDInt *diagslots, UCFDReal *values)
{
    Precon pc = (Precon)calloc(1, sizeof(*pc));
    UCFDCheckNull(pc, "Precon allocation failed\n");

    pc->type_name   = NULL;
    pc->rowptr      = rowptr;
    pc->colidx      = colidx;
    pc->diagslots   = diagslots;
    pc->values      = values;
    pc->data        = NULL;
    *precon         = pc;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPreconCreatefromMPIMat(Precon *precon, SpMat mat, UCFDInt n)
{
    Precon pc = (Precon)calloc(1, sizeof(*pc));
    UCFDCheckNull(pc, "Precon allocation failed\n");

    pc->type_name   = NULL;
    UCFDCall(UCFDMatCopyPattern(mat, &pc->rowptr, &pc->colidx, &pc->values));

    /* Construct local-pattern diagslots */
    pc->diagslots       = malloc((size_t)n * sizeof(UCFDInt));
    set_diagslots(n, pc->rowptr, pc->colidx, pc->diagslots);
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
    free(*precon);
    *precon = NULL;
    
    UCFDFunctionReturn(UCFD_SUCCESS);
}
