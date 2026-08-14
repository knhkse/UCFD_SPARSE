#include "ilu.h"
#include "inverse.h"


static ucfd_status_t PBILUPreconPreparePerColor(UCFDInt nstart, UCFDInt nend, Precon precon)
{
    Precon_BILU *bilu = (Precon_BILU *)precon->data;
    const UCFDInt block = bilu->block;
    const UCFDInt blkdim = block*block;
    const UCFDInt *restrict rowptr = precon->rowptr, *restrict colidx = precon->colidx, \
                  *restrict diagslots = precon->diagslots;
    UCFDInt *restrict iw = bilu->iw;
    UCFDReal *restrict values = precon->values;

    UCFDInt idx, kdx, row, col, ele;
    UCFDInt ck, kk, kst, ked, jj, iwj;
    UCFDReal v, Aik[block][block];

    OMPWrapper(kdx, row, col, ele, ck, kk, kst, ked, jj, iwj, v, Aik)
    for (idx=nstart; idx<nend; ++idx)
    {
        const UCFDInt st = rowptr[idx];
        const UCFDInt ed = diagslots[idx];
        const UCFDInt jed = rowptr[idx+1];

        for (kdx = st; kdx < ed; ++kdx)
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

static ucfd_status_t PBILUPreconPrepare(Precon precon)
{
    Precon_PBILU *pbilu = (Precon_PBILU *)precon->data;
    UCFDInt i;

    for (i=0; i<pbilu->ncolors; ++i)
        UCFDCall(PBILUPreconPreparePerColor(pbilu->icolors[i], pbilu->icolors[i+1], precon));
    
    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t PBILUPreconLowerApply(UCFDInt nstart, UCFDInt nend, Precon precon, UCFDReal *b)
{
    Precon_BILU *bilu = (Precon_BILU *)precon->data;
    const UCFDInt block = bilu->block;
    const UCFDInt blkdim = block*block;
    const UCFDInt *restrict rowptr = precon->rowptr, *restrict colidx = precon->colidx, \
                  *restrict diagslots = precon->diagslots;
    const UCFDReal *restrict values = precon->values;
    
    UCFDInt idx, jdx, kdx, row, col, cind;
    UCFDReal arr[block];

    // Forward substitution
    OMPWrapper(jdx, kdx, row, col, cind, arr)
    for (idx = nstart; idx < nend; ++idx)
    {
        const UCFDInt dd = diagslots[idx];
        const UCFDInt st = rowptr[idx];

        // Initialize arr
        for (kdx = 0; kdx < block; ++kdx)
            arr[kdx] = b[kdx + idx * block];

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

        for (kdx = 0; kdx < block; ++kdx)
            b[kdx + idx * block] = arr[kdx];
    }

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t PBILUPreconUpperApply(UCFDInt nstart, UCFDInt nend, Precon precon, UCFDReal *b)
{
    Precon_BILU *bilu = (Precon_BILU *)precon->data;
    const UCFDInt block = bilu->block;
    const UCFDInt blkdim = block*block;
    const UCFDInt *restrict rowptr = precon->rowptr, *restrict colidx = precon->colidx, \
                  *restrict diagslots = precon->diagslots;
    const UCFDReal *restrict values = precon->values;
    
    UCFDInt idx, jdx, kdx, row, col, cind;
    UCFDReal arr[block];

    // Backward substitution
    OMPWrapper(jdx, kdx, row, col, cind, arr)
    for (idx = nstart; idx < nend; ++idx)
    {
        const UCFDInt dd = diagslots[idx];
        const UCFDInt ed = rowptr[idx + 1];

        // Initialize
        for (kdx = 0; kdx < block; ++kdx)
            arr[kdx] = b[kdx + idx * block];

        for (jdx = dd + 1; jdx < ed; ++jdx)
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

        // LU substitution for vector
        lusub(block, &(values[dd*blkdim]), arr);
        for (row=0; row<block; ++row) b[idx*block+row] = arr[row];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t PBILUPreconApply(Precon precon, UCFDReal *b)
{
    Precon_PBILU *pbilu = (Precon_PBILU *)precon->data;
    UCFDInt i;

    // Lower sweep
    for (i=0; i<pbilu->ncolors; ++i)
        UCFDCall(PBILUPreconLowerApply(pbilu->icolors[i], pbilu->icolors[i+1], precon, b));
    
    // Upper sweep
    for (i=pbilu->ncolors-1; i>=0; --i)
        UCFDCall(PBILUPreconUpperApply(pbilu->icolors[i], pbilu->icolors[i+1], precon, b));
    
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPreconSetPBILU(Precon *precon, UCFDInt bn, UCFDInt block, UCFDInt ncolors, UCFDInt *icolors)
{
    UCFDCheckNull(*precon, "Preconditioner must be initialized\n");
    Precon pc = *precon;
    Precon_PBILU *pbilu = (Precon_PBILU *)calloc(1, sizeof(*pbilu));
    UCFDCheckNull(pbilu, "PBILU precon allocation failed\n");

    ((Precon_BILU *)pbilu)->bn          = bn;
    ((Precon_BILU *)pbilu)->block       = block;
    ((Precon_BILU *)pbilu)->iw          = (UCFDInt *)malloc((size_t)bn*sizeof(UCFDInt));
    pbilu->ncolors                      = ncolors;
    pbilu->icolors                      = icolors;

    /* Initialize working array */
    for (UCFDInt i = 0; i < bn; ++i) ((Precon_BILU *)pbilu)->iw[i] = -1;

    pc->type_name                       = PBILU;
    pc->data                            = pbilu;
    pc->ops->prepare                    = PBILUPreconPrepare;
    pc->ops->apply                      = PBILUPreconApply;
    pc->ops->destroy                    = ILUPreconDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
