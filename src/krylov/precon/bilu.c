#include "ilu.h"
#include "inverse.h"


static ucfd_status_t BILUPreconPrepare(Precon precon)
{
    Precon_BILU *bilu = (Precon_BILU *)precon->data;
    const UCFDInt bn = bilu->bn, block = bilu->block;
    const UCFDInt blkdim = block*block;
    const UCFDInt *restrict rowptr = precon->rowptr, *restrict colidx = precon->colidx, \
                  *restrict diagslots = precon->diagslots;
    UCFDInt *restrict iw = bilu->iw;
    UCFDReal *restrict values = precon->values;

    UCFDInt idx, kdx, row, col, ele;
    UCFDInt ck, kk, kst, ked, jj, iwj;
    UCFDReal v, Aik[block][block];

    for (idx=0; idx<bn; ++idx)
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
            
            // A[i,k] := A[i,k] @ inv(A[k,k])
            lusubmattrans(block, &(values[kk*blkdim]), &(values[kdx*blkdim]));
            // memcpy(Aik, &values[kdx*blkdim], sizeof(double)*block);
            for (row=0; row<block; ++row) {
                for (col=0; col<block; ++col)
                    Aik[row][col] = values[kdx*blkdim+row*block+col];
            }

            // Prepare iw
            for (jj=kst; jj<ked; ++jj) iw[colidx[jj]] = jj;

            // j iteration
            for (jj=kdx+1; jj<jed; ++jj) {
                iwj = iw[colidx[jj]];

                if (iwj != -1) {
                    // values[jj] -= Aik * values[iwj]
                    for (row=0; row<block; ++row) {
                        for (col=0; col<block; ++col) {
                            v = 0.0;
                            for (ele=0; ele<block; ++ele)
                                v += Aik[row][ele] * values[iwj*blkdim+ele*block+col];
                            values[jj*blkdim+row*block+col] -= v;
                        }
                    }
                }
            }

            // Clean iw
            for (jj=kst; jj<ked; ++jj) iw[colidx[jj]] = -1;
        }
        // LU decomposition of current row diagonal matrix
        ludcmp(block, &(values[ed*blkdim]));
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t BILUPreconApply(Precon precon, UCFDReal *b)
{
    Precon_BILU *bilu = (Precon_BILU *)precon->data;
    const UCFDInt bn = bilu->bn, block = bilu->block;
    const UCFDInt blkdim = block*block;
    const UCFDInt *restrict rowptr = precon->rowptr, *restrict colidx = precon->colidx, \
                  *restrict diagslots = precon->diagslots;
    const UCFDReal *restrict values = precon->values;

    UCFDInt idx, jdx, kdx, row, col, cind;
    UCFDReal arr[block];

    // Forward substitution
    for (idx=0; idx<bn; ++idx)
    {
        const UCFDInt dd = diagslots[idx];
        const UCFDInt st = rowptr[idx];

        // Initialize arr
        for (kdx=0; kdx<block; ++kdx)
            arr[kdx] = b[kdx + idx * block];

        for (jdx=st; jdx<dd; ++jdx)
        {
            cind = colidx[jdx];

            for (row=0; row<block; ++row)
            {
                UCFDReal v = 0.0;
                for (col = 0; col < block; ++col)
                    v += values[col + row * block + jdx * blkdim] * b[col + cind * block];
                arr[row] -= v;
            }
        }

        for (kdx=0; kdx<block; ++kdx)
            b[kdx + idx * block] = arr[kdx];
    }

    // Backward substitution
    for (idx=bn-1; idx>-1; --idx)
    {
        // Initialize
        for (kdx = 0; kdx < block; ++kdx)
            arr[kdx] = b[kdx + idx * block];

        const UCFDInt dd = diagslots[idx];
        const UCFDInt ed = rowptr[idx + 1];

        for (jdx=dd+1; jdx<ed; ++jdx)
        {
            cind = colidx[jdx];

            for (row=0; row<block; ++row)
            {
                UCFDReal v = 0.0;
                for (col=0; col<block; ++col)
                    v += values[col + row * block + jdx * blkdim] * b[col + cind * block];
                arr[row] -= v;
            }
        }

        // LU substitution for vector
        lusub(block, &(values[dd*blkdim]), arr);
        for (row=0; row<block; ++row) b[idx*block+row] = arr[row];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPreconSetBILU(Precon *precon, UCFDInt bn, UCFDInt block)
{
    UCFDCheckNull(*precon, "Preconditioner must be initialized\n");
    Precon pc = *precon;
    Precon_BILU *bilu = (Precon_BILU *)calloc(1, sizeof(*bilu));
    UCFDCheckNull(bilu, "BILU precon allocation failed\n");

    bilu->bn            = bn;
    bilu->block         = block;
    bilu->iw            = (UCFDInt *)malloc((size_t)bn*sizeof(UCFDInt));
    
    /* Initialize working array */
    for (UCFDInt i = 0; i < bn; ++i) bilu->iw[i] = -1;

    pc->type_name       = BILU;
    pc->data            = bilu;
    pc->ops->prepare    = BILUPreconPrepare;
    pc->ops->apply      = BILUPreconApply;
    pc->ops->destroy    = ILUPreconDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}