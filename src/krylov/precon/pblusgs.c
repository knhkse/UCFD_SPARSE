#include <string.h>

#include "blusgs.h"
#include "inverse.h"


static ucfd_status_t PBLUSGSPreconLowerApply(UCFDInt nstart, UCFDInt nend, Precon precon, UCFDReal *b)
{
    Precon_BLUSGS *blu = (Precon_BLUSGS *)precon->data;
    const UCFDInt bn = blu->bn, block = blu->block;
    const UCFDInt blkdim = block*block;
    const UCFDInt *restrict rowptr = precon->rowptr, *restrict colidx = precon->colidx, \
                  *restrict diagslots = precon->diagslots;
    const UCFDReal *restrict values = precon->values, *restrict diagvalues = blu->diagvalues;

    UCFDInt idx, jdx, kdx, row, col;
    UCFDInt dd, st, cind;
    UCFDReal v, arr[block];

    // Forward substitution
    OMPWrapper(jdx, kdx, row, col, dd, st, cind, v, arr)
    for (idx=nstart; idx<nend; ++idx)
    {
        st = rowptr[idx];
        dd = diagslots[idx];

        for (kdx=0; kdx<block; ++kdx)
            arr[kdx] = b[kdx + idx * block];

        for (jdx=st; jdx<dd; ++jdx)
        {
            cind = colidx[jdx];
            for (row=0; row<block; ++row) {
                v = 0.0;
                for (col=0; col<block; ++col)
                    v += values[col+row*block+jdx*blkdim] * b[col+cind*block];
                arr[row] -= v;
            }
        }

        lusub(block, &diagvalues[idx+blkdim], arr);
        for (kdx=0; kdx<block; ++kdx)
            b[kdx + idx*block] = arr[kdx];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t PBLUSGSPreconUpperApply(UCFDInt nstart, UCFDInt nend, Precon precon, UCFDReal *b)
{
    Precon_BLUSGS *blu = (Precon_BLUSGS *)precon->data;
    const UCFDInt bn = blu->bn, block = blu->block;
    const UCFDInt blkdim = block*block;
    const UCFDInt *restrict rowptr = precon->rowptr, *restrict colidx = precon->colidx, \
                  *restrict diagslots = precon->diagslots;
    const UCFDReal *restrict values = precon->values, *restrict diagvalues = blu->diagvalues;

    UCFDInt idx, jdx, kdx, row, col;
    UCFDInt dd, ed, cind;
    UCFDReal v, arr[block];

    // Backward substitution
    OMPWrapper(jdx, kdx, row, col, dd, ed, cind, v, arr)
    for (idx=nstart; idx<nend; ++idx)
    {
        dd = diagslots[idx];
        ed = rowptr[idx+1];

        for (kdx = 0; kdx < block; ++kdx)
            arr[kdx] = 0.0;

        for (jdx=dd+1; jdx<ed; ++jdx)
        {
            cind = colidx[jdx];
            for (row=0; row<block; ++row) {
                v = 0.0;
                for (col=0; col<block; ++col)
                    v += values[col+row*block+jdx*blkdim]*b[col+cind*block];
                arr[row] += v;
            }
        }
        lusub(block, &diagvalues[idx*blkdim], arr);

        for (kdx=0; kdx<block; ++kdx)
            b[kdx + idx*block] -= arr[kdx];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t PBLUSGSPreconApply(Precon precon, UCFDReal *b)
{
    Precon_PBLUSGS *blu = (Precon_PBLUSGS *)precon->data;
    UCFDInt i;
    const UCFDInt ncolors = blu->ncolors;
    const UCFDInt *icolors = blu->icolors;

    for (i=0; i<ncolors; ++i)
        UCFDCall(PBLUSGSPreconLowerApply(icolors[i], icolors[i+1], precon, b));
    
    for (i=ncolors-1; i>=0; --i)
        UCFDCall(PBLUSGSPreconUpperApply(icolors[i], icolors[i+1], precon, b));

    UCFDFunctionReturn(UCFD_SUCCESS);
}


ucfd_status_t UCFDPreconSetPBLUSGS(Precon *precon, UCFDInt bn, UCFDInt block, UCFDInt ncolors, UCFDInt *icolors)
{
    UCFDCheckNull(*precon, "Preconditioner must be initialized\n");
    Precon pc = *precon;
    Precon_PBLUSGS *pblu = (Precon_PBLUSGS *)calloc(1, sizeof(*pblu));
    UCFDCheckNull(pblu, "PBLUSGS precon allocation failed\n");

    ((Precon_BLUSGS *)pblu)->bn             = bn;
    ((Precon_BLUSGS *)pblu)->block          = block;
    ((Precon_BLUSGS *)pblu)->diagvalues     = (UCFDReal *)calloc((size_t)bn*block*block, sizeof(UCFDReal));
    pblu->ncolors                           = ncolors;
    pblu->icolors                           = icolors;

    pc->type_name       = PBLUSGS;
    pc->data            = pblu;
    pc->ops->prepare    = BLUSGSPreconPrepare;
    pc->ops->apply      = PBLUSGSPreconApply;
    pc->ops->destroy    = BLUSGSPreconDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
