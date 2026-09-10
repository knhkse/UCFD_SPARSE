#include "ilu.h"


static void _ilu_prepare(const UCFDInt n,
                         const UCFDInt *restrict rowptr, const UCFDInt *restrict colidx,
                         const UCFDInt *restrict diagslots, UCFDInt *restrict iw,
                         UCFDReal *restrict values)
{
    UCFDInt idx, kdx, ck;
    UCFDInt kk, kst, ked, jj, iwj;
    UCFDReal Aik;

    for (idx=0; idx<n; ++idx)
    {
        const UCFDInt st = rowptr[idx];
        const UCFDInt ed = diagslots[idx];
        const UCFDInt jed = rowptr[idx+1];

        for (kdx=st; kdx<ed; ++kdx)
        {
            ck = colidx[kdx];
            kk = diagslots[ck];
            kst = rowptr[ck];
            ked = rowptr[ck+1];
            
            /* A[i,k] := A[i,k] / A[k,k] */
            values[kdx] /= values[kk];
            Aik = values[kdx];

            /* Prepare iw */
            for (jj = kst; jj < ked; ++jj) iw[colidx[jj]] = jj;

            /* j iteration */
            /* A[i,j] -= A[i,k]*A[k,k] */
            for (jj=kdx+1; jj<jed; ++jj) {
                iwj = iw[colidx[jj]];
                if (iwj != -1) values[jj] -= Aik*values[iwj];
            }

            /* Clean iw */
            for (jj=kst; jj<ked; ++jj) iw[colidx[jj]] = -1;
        }
    }
}

static void _ilu_apply(const UCFDInt n,
                       const UCFDInt *restrict rowptr, const UCFDInt *restrict colidx,
                       const UCFDInt *restrict diagslots, UCFDReal *restrict values,
                       UCFDReal *restrict b)
{
    UCFDInt idx, jdx, cind;

    /* Forward sweep */
    for (idx=0; idx<n; ++idx)
    {
        const UCFDInt dd = diagslots[idx];
        const UCFDInt st = rowptr[idx];
        UCFDReal v = b[idx];

        for (jdx=st; jdx<dd; ++jdx)
        {
            cind = colidx[jdx];
            v -= values[jdx] * b[cind];
        }
        b[idx] = v;
    }

    /* Backward sweep */
    for (idx=n-1; idx>-1; --idx)
    {
        const UCFDInt dd = diagslots[idx];
        const UCFDInt ed = rowptr[idx+1];
        UCFDReal v = b[idx];

        for (jdx=dd; jdx<ed; ++jdx)
        {
            cind = colidx[jdx];
            v -= values[jdx] * b[cind];
        }
        b[idx] = v/values[dd];
    }
}


static ucfd_status_t ILUPreconPrepare(Precon precon)
{
    Precon_ILU *ilu = (Precon_ILU *)precon->data;
    _ilu_prepare(
        ilu->n, precon->rowptr, precon->colidx, precon->diagslots,
        ilu->iw, precon->values
    );
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t ILUPreconApply(Precon precon, UCFDReal *b)
{
    Precon_ILU *ilu = (Precon_ILU *)precon->data;
    _ilu_apply(
        ilu->n, precon->rowptr, precon->colidx, precon->diagslots,
        precon->values, b
    );
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
    for (UCFDInt i=0; i<n; ++i) ilu->iw[i] = -1;

    pc->type_name       = ILU;
    pc->values          = malloc(pc->nnz*sizeof(UCFDReal));
    pc->data            = ilu;
    pc->ops->prepare    = ILUPreconPrepare;
    pc->ops->apply      = ILUPreconApply;
    pc->ops->destroy    = ILUPreconDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}