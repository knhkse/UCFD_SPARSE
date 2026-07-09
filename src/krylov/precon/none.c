#include "ucfdpcimpl.h"


static ucfd_status_t NonePreconFunction(Precon precon)
{
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t NonePreconApply(Precon precon, UCFDReal *x)
{
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDCreateNonePrecon(Precon *precon)
{
    Precon pc = (Precon)calloc(1, sizeof(*pc));
    UCFDNullCheck(pc, "Precon allocation failed\n");

    pc->type_name       = NONE;
    pc->rowptr          = NULL;
    pc->colidx          = NULL;
    pc->diagslots       = NULL;
    pc->data            = NULL;
    pc->ops->prepare    = NonePreconFunction;
    pc->ops->apply      = NonePreconApply;
    pc->ops->destroy    = NonePreconFunction;

    *precon = pc;

    UCFDFunctionReturn(UCFD_SUCCESS);
}