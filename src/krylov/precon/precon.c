/**
 * General functions for preconditioner
 */
#include <stdio.h>
#include "ucfdpcimpl.h"
#include "ucfdmatimpl.h"


char preconwarning[] =  "Preconditioner warning :\n"
                        "Same pointer address is indicated in both system matrix and preconditioner matrix.\n"
                        "If ILU-type preconditioner is used, incorrect result can be occurred.\n";

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

ucfd_status_t UCFDPreconPrepare(Precon precon)
{
    UCFDCheckNull(precon->type_name, "Preconditioner type must be set\n");
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
