#include <string.h>

#include "blusgs.h"
#include "inverse.h"


ucfd_status_t BLUSGSPreconPrepare(Precon precon)
{
    Precon_BLUSGS *blu = (Precon_BLUSGS *)precon->data;
    const UCFDInt *restrict diagslots = precon->diagslots;
    const UCFDReal *restrict values = precon->values;
    UCFDReal *restrict diagvalues = blu->diagvalues;
    const UCFDInt bn = blu->bn, block = blu->block;
    const UCFDInt blkdim = block*block;
    UCFDInt idx, didx;
    UCFDReal *diagblock;

    OMPWrapper(didx, diagblock)
    for (idx=0; idx<bn; ++idx) {
        didx = diagslots[idx];
        diagblock = &diagvalues[idx*blkdim];
        memcpy(diagblock, &values[didx*blkdim], sizeof(UCFDReal)*blkdim);
        ludcmp(block, diagblock);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t BLUSGSPreconApply(Precon precon, UCFDReal *b)
{
    Precon_BLUSGS *blu = (Precon_BLUSGS *)precon->data;
    const UCFDInt *restrict rowptr = precon->rowptr, *restrict colidx = precon->colidx, \
                  *restrict diagslots = precon->diagslots;
    const UCFDReal *restrict values = precon->values, *restrict diagvalues = blu->diagvalues;
    const UCFDInt bn = blu->bn, block = blu->block;
    const UCFDInt blkdim = block*block;

    UCFDInt idx, jdx, kdx, row, col, cind;
    UCFDReal arr[block];

    // Forward sweep : (D+L)x' = b -> x' = inv(D) * (b-Lx')
    for (idx = 0; idx < bn; ++idx)
    {
        const UCFDInt dd = diagslots[idx];
        const UCFDInt st = rowptr[idx];

        // arr := b
        for (kdx = 0; kdx < block; ++kdx)
            arr[kdx] = b[kdx + idx * block];

        // arr := b - Lx'
        for (jdx = st; jdx < dd; ++jdx)
        {
            cind = colidx[jdx];
            for (row = 0; row < block; ++row)
            {
                UCFDReal v = 0.0;
                for (col = 0; col < block; ++col)
                    v += values[col + row * block + jdx * blkdim] * b[col + cind * block];
                arr[row] -= v;
            }
        }

        // x' := inv(D) * (b-Lx') = inv(D) * arr
        lusub(block, &diagvalues[idx * blkdim], arr);
        for (kdx = 0; kdx < block; ++kdx)
            b[kdx + idx * block] = arr[kdx];
    }

    // Backward sweep : (D+U)x = Dx' -> x = x' - inv(D) * Ux
    for (idx = bn - 1; idx > -1; --idx)
    {
        const UCFDInt dd = diagslots[idx];
        const UCFDInt ed = rowptr[idx + 1];

        // Initialize
        for (kdx = 0; kdx < block; ++kdx)
            arr[kdx] = 0.0;

        // arr := Ux
        for (jdx = dd + 1; jdx < ed; ++jdx)
        {
            cind = colidx[jdx];
            for (row = 0; row < block; ++row)
            {
                UCFDReal v = 0.0;
                for (col = 0; col < block; ++col)
                    v += values[col + row * block + jdx * blkdim] * b[col + cind * block];
                arr[row] += v;
            }
        }

        // arr := inv(D) Ux
        lusub(block, &diagvalues[idx * blkdim], arr);

        // b := b - inv(D) Ux
        for (kdx = 0; kdx < block; ++kdx)
            b[kdx + idx * block] -= arr[kdx];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t BLUSGSPreconDestroy(Precon precon)
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