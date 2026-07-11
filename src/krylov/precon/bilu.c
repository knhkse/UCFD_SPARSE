#include "bilu.h"
#include "inverse.h"


static ucfd_status_t BILUPreconPrepare(Precon precon)
{
    Precon_BILU *bilu = (Precon_BILU *)precon->data;
    UCFDInt bn = bilu->bn;
    UCFDInt block = bilu->block;
    UCFDInt idx, kdx, ck, row, col, ele;
    UCFDInt st, ed, jed, kk, kst, ked, jj, iwj;
    UCFDReal v, Aik[block][block];
    UCFDInt blkdim = block*block;

    for (idx=0; idx<bn; idx++)
    {
        st = precon->rowptr[idx];
        ed = precon->diagslots[idx];
        jed = precon->rowptr[idx+1];

        for (kdx = st; kdx < ed; kdx++)
        {
            ck = precon -> colidx[kdx];
            kk = precon -> diagslots[ck];
            kst = precon -> rowptr[ck];
            ked = precon -> rowptr[ck+1];
            
            // A[i,k] := A[i,k] @ inv(A[k,k])
            lusubmattrans(block, &(precon->values[kk*blkdim]), &(precon->values[kdx*blkdim]));
            // memcpy(Aik, &values[kdx*blkdim], sizeof(double)*block);
            for (row=0; row<block; row++) {
                for (col=0; col<block; col++)
                    Aik[row][col] = precon->values[kdx*blkdim+row*block+col];
            }

            // Prepare iw
            for (jj=kst; jj<ked; jj++) bilu->iw[precon -> colidx[jj]] = jj;

            // j iteration
            for (jj=kdx+1; jj<jed; jj++) {
                iwj = bilu->iw[precon -> colidx[jj]];

                if (iwj != -1) {
                    // values[jj] -= Aik * values[iwj]
                    for (row=0; row<block; row++) {
                        for (col=0; col<block; col++) {
                            v = 0.0;
                            for (ele=0; ele<block; ele++)
                                v += Aik[row][ele] * precon->values[iwj*blkdim+ele*block+col];
                            precon->values[jj*blkdim+row*block+col] -= v;
                        }
                    }
                }
            }

            // Clean iw
            for (jj=kst; jj<ked; jj++) bilu->iw[precon -> colidx[jj]] = -1;
        }
        // LU decomposition of current row diagonal matrix
        ludcmp(block, &(precon->values[ed*blkdim]));
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t BILUPreconApply(Precon precon, UCFDReal *b)
{
    Precon_BILU *bilu = (Precon_BILU *)precon->data;
    UCFDInt bn = bilu->bn;
    UCFDInt block = bilu->block;
    UCFDInt idx, jdx, kdx, row, col;
    UCFDInt dd, st, ed, cind;
    UCFDReal v, arr[block];
    UCFDInt blkdim = block*block;

    // Forward substitution
    for (idx = 0; idx < bn; idx++)
    {
        dd = precon->diagslots[idx];
        st = precon->rowptr[idx];

        // Initialize arr
        for (kdx = 0; kdx < block; kdx++)
            arr[kdx] = b[kdx + idx * block];

        for (jdx = st; jdx < dd; jdx++)
        {
            cind = precon -> colidx[jdx];

            for (row = 0; row < block; row++)
            {
                v = 0.0;
                for (col = 0; col < block; col++)
                    v += precon->values[col + row * block + jdx * blkdim] * b[col + cind * block];
                arr[row] -= v;
            }
        }

        for (kdx = 0; kdx < block; kdx++)
            b[kdx + idx * block] = arr[kdx];
    }

    // Backward substitution
    for (idx = bn - 1; idx > -1; idx--)
    {
        dd = precon->diagslots[idx];
        ed = precon->rowptr[idx + 1];

        // Initialize
        for (kdx = 0; kdx < block; kdx++)
            arr[kdx] = b[kdx + idx * block];

        for (jdx = dd + 1; jdx < ed; jdx++)
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

        // LU substitution for vector
        lusub(block, &(precon->values[dd*blkdim]), arr);
        for (row=0; row<block; row++) b[idx*block+row] = arr[row];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t BILUPreconDestroy(Precon precon)
{
    if (!precon) UCFDFunctionReturn(UCFD_SUCCESS);
    Precon_BILU *bilu = (Precon_BILU *)precon->data;
    free(bilu->iw);
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPreconSetBILU(Precon *precon, UCFDInt bn, UCFDInt block)
{
    Precon pc = *precon;
    Precon_BILU *bilu = (Precon_BILU *)calloc(1, sizeof(*bilu));
    UCFDCheckNull(bilu, "BILU precon allocation failed\n");

    bilu->bn            = bn;
    bilu->block         = block;
    bilu->iw            = (UCFDInt *)malloc((size_t)bn*sizeof(UCFDInt));
    
    /* Initialize working array */
    for (UCFDInt i = 0; i < bn; i++) bilu->iw[i] = -1;

    pc->type_name       = BILU;
    pc->data            = bilu;
    pc->ops->prepare    = BILUPreconPrepare;
    pc->ops->apply      = BILUPreconApply;
    pc->ops->destroy    = BILUPreconDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}