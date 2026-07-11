#include <string.h>

#include "blusgs.h"
#include "inverse.h"


static ucfd_status_t BLUSGSPreconPrepare(Precon precon)
{
    Precon_BLUSGS *blu = (Precon_BLUSGS *)precon->data;
    UCFDInt bn = blu->bn;
    UCFDInt block = blu->block;
    UCFDInt idx, didx;
    UCFDInt blkdim = block*block;
    UCFDReal *diagblock;

    OMPWrapper(didx, diagblock)
    for (idx=0; idx<bn; idx++) {
        didx = precon->diagslots[idx];
        diagblock = &blu->diagvalues[idx*blkdim];
        memcpy(diagblock, &precon->values[didx*blkdim], sizeof(UCFDReal)*blkdim);
        ludcmp(block, diagblock);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t BLUSGSPreconApply(Precon precon, UCFDReal *b)
{
    Precon_BLUSGS *blu = (Precon_BLUSGS *)precon->data;
    UCFDInt bn = blu->bn;
    UCFDInt block = blu->block;
    UCFDInt idx, jdx, kdx, row, col;
    UCFDInt dd, st, ed, cind;
    UCFDReal v, arr[block];
    UCFDInt blkdim = block*block;

    // Forward sweep : (D+L)x' = b -> x' = inv(D) * (b-Lx')
    for (idx = 0; idx < bn; idx++)
    {
        dd = precon->diagslots[idx];
        st = precon->rowptr[idx];

        // arr := b
        for (kdx = 0; kdx < block; kdx++)
            arr[kdx] = b[kdx + idx * block];

        // arr := b - Lx'
        for (jdx = st; jdx < dd; jdx++)
        {
            cind = precon->colidx[jdx];
            for (row = 0; row < block; row++)
            {
                v = 0.0;
                for (col = 0; col < block; col++)
                    v += precon->values[col + row * block + jdx * blkdim] * b[col + cind * block];
                arr[row] -= v;
            }
        }

        // x' := inv(D) * (b-Lx') = inv(D) * arr
        lusub(block, &blu->diagvalues[idx * blkdim], arr);
        for (kdx = 0; kdx < block; kdx++)
            b[kdx + idx * block] = arr[kdx];
    }

    // Backward sweep : (D+U)x = Dx' -> x = x' - inv(D) * Ux
    for (idx = bn - 1; idx > -1; idx--)
    {
        dd = precon->diagslots[idx];
        ed = precon->rowptr[idx + 1];

        // Initialize
        for (kdx = 0; kdx < block; kdx++)
            arr[kdx] = 0.0;

        // arr := Ux
        for (jdx = dd + 1; jdx < ed; jdx++)
        {
            cind = precon->colidx[jdx];
            for (row = 0; row < block; row++)
            {
                v = 0.0;
                for (col = 0; col < block; col++)
                    v += precon->values[col + row * block + jdx * blkdim] * b[col + cind * block];
                arr[row] += v;
            }
        }

        // arr := inv(D) Ux
        lusub(block, &blu->diagvalues[idx * blkdim], arr);

        // b := b - inv(D) Ux
        for (kdx = 0; kdx < block; kdx++)
            b[kdx + idx * block] -= arr[kdx];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t BLUSGSPreconDestroy(Precon precon)
{
    if (!precon) UCFDFunctionReturn(UCFD_SUCCESS);
    Precon_BLUSGS *blu = (Precon_BLUSGS *)precon->data;
    free(blu->diagvalues);
    UCFDFunctionReturn(UCFD_SUCCESS);
}


ucfd_status_t UCFDPreconSetBLUSGS(Precon *precon, UCFDInt bn, UCFDInt block)
{
    UCFDCheckNull(*precon, "Preconditioner must be initialized\n");
    Precon pc = *precon;
    Precon_BLUSGS *blu = (Precon_BLUSGS *)calloc(1, sizeof(*blu));
    UCFDCheckNull(blu, "BLU-SGS precon allocation failed\n");

    blu->bn             = bn;
    blu->block          = block;
    blu->diagvalues     = (UCFDReal *)calloc((size_t)bn*block*block, sizeof(UCFDReal));

    pc->type_name       = BLUSGS;
    pc->data            = blu;
    pc->ops->prepare    = BLUSGSPreconPrepare;
    pc->ops->apply      = BLUSGSPreconApply;
    pc->ops->destroy    = BLUSGSPreconDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}