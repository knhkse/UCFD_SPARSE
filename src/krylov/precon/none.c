#include "ucfdpcimpl.h"


static inline ucfd_status_t NonePreconFunction(Precon precon)
{
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static inline ucfd_status_t NonePreconApply(Precon precon, UCFDReal *x)
{
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPreconSetNone(Precon *precon)
{
    UCFDCheckNull(*precon, "Preconditioner must be initialized\n");
    Precon pc = *precon;
    UCFDCheckNull(pc, "Precon allocation failed\n");

    pc->type_name       = NONE;
    pc->values          = NULL;
    pc->data            = NULL;
    pc->ops->prepare    = NonePreconFunction;
    pc->ops->apply      = NonePreconApply;
    pc->ops->destroy    = NonePreconFunction;

    *precon = pc;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
