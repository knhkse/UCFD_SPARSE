#include "ucfdpcimpl.h"


static inline ucfd_status_t NonePreconFunction(Precon precon)
{
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline ucfd_status_t NonePreconApply(Precon precon, UCFDReal *x)
{
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPreconCreateNone(Precon *precon)
{
    Precon pc = (Precon)calloc(1, sizeof(*pc));
    UCFDCheckNull(pc, "Precon allocation failed\n");

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