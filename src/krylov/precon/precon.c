/**
 * General functions for preconditioner
 */

#include "ucfdpcimpl.h"
#include "ucfdmatimpl.h"


ucfd_status_t UCFDPreconInitfromMatrix(Precon *precon, SpMat sysmat, UCFDInt *diagslots)
{
    UCFDNullCheck(sysmat->type_name, "Input matrix is empty\n");
    Precon pc = (Precon)calloc(1, sizeof(*pc));
    UCFDNullCheck(pc, "Precon allocation failed\n");

    pc->type_name   = NULL;
    pc->rowptr      = sysmat->rowptr;
    pc->colidx      = sysmat->colidx;
    pc->diagslots   = diagslots;
    pc->data        = NULL;
    *precon         = pc;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPreconInitfromArrays(Precon *precon, UCFDInt *rowptr, UCFDInt *colidx, UCFDInt *diagslots)
{
    Precon pc = (Precon)calloc(1, sizeof(*pc));
    UCFDNullCheck(pc, "Precon allocation failed\n");

    pc->type_name   = NULL;
    pc->rowptr      = rowptr;
    pc->colidx      = colidx;
    pc->diagslots   = diagslots;
    pc->data        = NULL;
    *precon         = pc;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPreconPrep(Precon precon)
{
    UCFDNullCheck(precon->type_name, "Preconditioner type must be set\n");
    UCFDCall(precon->ops->prepare(precon));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


ucfd_status_t UCFDPreconDestroy(Precon *precon)
{
    if (!precon || !*precon) UCFDFunctionReturn(UCFD_SUCCESS);
    UCFDCall((*precon)->ops->destroy(*precon));
    free((*precon)->data);
    free(*precon);
    *precon = NULL;
    
    UCFDFunctionReturn(UCFD_SUCCESS);
}