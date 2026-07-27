#include "ilu.h"


static ucfd_status_t ILUPreconPrepare(Precon precon)
{
    Precon_ILU *ilu = (Precon_ILU *)precon->data;
    UCFDInt n = ilu->n;
    UCFDInt idx, kdx, ck;
    UCFDInt st, ed, jed, kk, kst, ked, jj, iwj;
    UCFDReal Aik;

    for (idx=0; idx<n; idx++)
    {
        st = precon->rowptr[idx];
        ed = precon->diagslots[idx];
        jed = precon->rowptr[idx+1];

        for (kdx=st; kdx<ed; kdx++)
        {
            ck = precon -> colidx[kdx];
            kk = precon -> diagslots[ck];
            kst = precon -> rowptr[ck];
            ked = precon -> rowptr[ck+1];
            
            /* A[i,k] := A[i,k] / A[k,k] */
            precon->values[kdx] /= precon->values[kk];
            Aik = precon->values[kdx];

            /* Prepare iw */
            for (jj = kst; jj < ked; jj++) ilu->iw[precon->colidx[jj]] = jj;

            /* j iteration */
            /* A[i,j] -= A[i,k]*A[k,k] */
            for (jj=kdx+1; jj<jed; jj++) {
                iwj = ilu->iw[precon->colidx[jj]];
                if (iwj != -1) precon->values[jj] -= Aik*precon->values[iwj];
            }

            /* Clean iw */
            for (jj=kst; jj<ked; jj++) ilu->iw[precon->colidx[jj]] = -1;
        }
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t ILUPreconApply(Precon precon, UCFDReal *b)
{
    Precon_ILU *ilu = (Precon_ILU *)precon->data;
    UCFDInt n = ilu->n;
    UCFDInt idx, jdx, dd, st, ed, cind;
    UCFDReal v;

    /* Forward sweep */
    for (idx=0; idx<n; idx++)
    {
        dd = precon->diagslots[idx];
        st = precon->rowptr[idx];

        v = b[idx];
        for (jdx=st; jdx<dd; jdx++)
        {
            cind = precon->colidx[jdx];
            v -= precon->values[jdx] * b[cind];
        }
        b[idx] = v;
    }

    /* Backward sweep */
    for (idx=n-1; idx>-1; idx--)
    {
        dd = precon->diagslots[idx];
        ed = precon->rowptr[idx+1];
        v = b[idx];

        for (jdx=dd; jdx<ed; jdx++)
        {
            cind = precon->colidx[jdx];
            v -= precon->values[jdx] * b[cind];
        }
        b[idx] = v/precon->values[dd];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t ILUPreconDestroy(Precon precon)
{
    if (!precon) UCFDFunctionReturn(UCFD_SUCCESS);
    Precon_ILU *ilu = (Precon_ILU *)precon->data;
    free(ilu->iw);
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPreconSetILU(Precon *precon, UCFDInt n)
{
    Precon pc = *precon;
    Precon_ILU *ilu = (Precon_ILU *)calloc(1, sizeof(*ilu));
    UCFDCheckNull(ilu, "ILU precon allocation failed\n");

    ilu->n      = n;
    ilu->iw     = (UCFDInt *)malloc((size_t)n*sizeof(UCFDInt));

    /* Initialize */
    for (UCFDInt i=0; i<n; i++) ilu->iw[i] = -1;

    pc->type_name       = ILU;
    pc->data            = ilu;
    pc->ops->prepare    = ILUPreconPrepare;
    pc->ops->apply      = ILUPreconApply;
    pc->ops->destroy    = ILUPreconDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}