#include <string.h>
#include "flowsys.h"
#include "flux.h"
#include "inverse.h"


/**
 * ! ----------------- array structure -------------
 * rhs, du, dup : [nlocal, nvars]   <- different with pyBaram
 * diag : [nlocal, nfvars, nfvars]
 * tdiag : [nlocal, nturbvars, nturbvars]
 */

/**
 * Pack & update kernels => per-element execution
 */
static ucfd_status_t
rank_blusgs_pack(const UCFDInt nlocal, const UCFDInt neles, const UCFDInt nvars,
                 const UCFDInt nfvars, const UCFDReal a0,
                 const UCFDInt *restrict cell_ids,
                 const UCFDReal *restrict rhs,
                 const UCFDReal *restrict dt,
                 UCFDReal *restrict rank_rhs,
                 UCFDReal *restrict diag)
{
    UCFDInt idx, ridx, kdx, row, col;
    const UCFDInt dim2 = nfvars*nfvars;

    for (idx=0; idx<neles; ++idx)
    {
        ridx = cell_ids[idx];
        for (kdx=0; kdx<nvars; ++kdx)
            // rank_rhs[ridx, kdx] = rhs[kdx, idx]
            rank_rhs[kdx + ridx*nvars] = rhs[idx + kdx*neles];
        
        for (row=0; row<nfvars; ++row) {
            const UCFDInt offset = row*nfvars + ridx*dim2;
            for (col=0; col<nfvars; ++col)
                // diag[ridx, row, col] = 0.0
                diag[col + offset] = 0.0;
            diag[row + offset] = 1/dt[idx] + a0;
        }
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_Pack(FlowSys sys, UCFDInt eidx, UCFDReal a0)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;
    FlowElem *e  = &sys->eles[eidx];
    
    UCFDCall(rank_blusgs_pack(
        sys->nlocal, e->neles, sys->nvars, sys->nfvars, a0,
        e->cell_ids, e->rhs, e->dt, blusgs->rhs, blusgs->diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t
rank_tblusgs_pack(const UCFDInt nlocal, const UCFDInt neles, const UCFDInt nvars,
                        const UCFDInt nfvars, const UCFDReal a0, const UCFDReal factor,
                        srcjacobian dsrcf,
                        const UCFDInt *restrict cell_ids,
                        const UCFDReal *restrict uptsb,
                        const UCFDReal *restrict dsrc,
                        const UCFDReal *restrict dt,
                        UCFDReal *restrict diag)
{
    UCFDInt idx, ridx, kdx, row, col;
    const UCFDInt nturbvars = nvars - nfvars;
    const UCFDInt dim2 = nturbvars*nturbvars;
    UCFDReal u[nvars], d[nvars];

    for (idx=0; idx<neles; ++idx)
    {
        ridx = cell_ids[idx];
        for (row=0; row<nturbvars; ++row) {
            for (col=0; col<nturbvars; ++col)
                diag[col + row*nturbvars + ridx*dim2] = 0.0;
        }

        // Prepare cell uf and dsrc
        for (kdx=0; kdx<nvars; ++kdx) {
            u[kdx] = uptsb[idx + kdx*neles];
            d[kdx] = dsrc[idx + kdx*neles];
        }
        dsrcf(nvars, nturbvars, u, &diag[ridx*dim2], d);

        for (row=0; row<nturbvars; ++row)
            diag[row + row*nturbvars + ridx*dim2] += 1/(dt[idx]*factor) + a0;
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_KWSST_Pack(FlowSys sys, UCFDInt eidx, UCFDReal turb_factor, UCFDReal a0)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;
    srcjacobian dsrcf = kwsst_src_jacobian;
    FlowElem *e  = &sys->eles[eidx];

    UCFDCall(rank_tblusgs_pack(
        sys->nlocal, e->neles, sys->nvars, sys->nfvars, a0, turb_factor,
        dsrcf, e->cell_ids, e->uptsb, e->dsrc, e->dt, blusgs->tdiag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_SA_Pack(FlowSys sys, UCFDInt eidx, UCFDReal turb_factor, UCFDReal a0)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;
    srcjacobian dsrcf = sa_src_jacobian;
    FlowElem *e  = &sys->eles[eidx];

    UCFDCall(rank_tblusgs_pack(
        sys->nlocal, e->neles, sys->nvars, sys->nfvars, a0, turb_factor,
        dsrcf, e->cell_ids, e->uptsb, e->dsrc, e->dt, blusgs->tdiag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


/**
 * Update
 */
static ucfd_status_t
rank_blusgs_update(const UCFDInt nlocal, const UCFDInt neles, const UCFDInt nvars,
                   const UCFDInt *restrict cell_ids,
                   const UCFDReal *restrict rank_du,
                   UCFDReal *restrict upts)
{
    UCFDInt idx, ridx, kdx;

    OMPWrapper(ridx, kdx)
    for (idx=0; idx<neles; ++idx)
    {
        ridx = cell_ids[idx];
        for (kdx=0; kdx<nvars; ++kdx)
            // upts[kdx, idx] += du[ridx, kdx]
            upts[idx + kdx*neles] += rank_du[kdx + ridx*nvars];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_Update(FlowSys sys, UCFDInt eidx)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;
    FlowElem *e  = &sys->eles[eidx];

    UCFDCall(rank_blusgs_update(
        sys->nlocal, e->neles, sys->nvars,
        e->cell_ids, blusgs->du, e->uptsb
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


/**
 * Inner iteration residual
 */
static ucfd_status_t
rank_sub_residual(const UCFDInt nlocal, const UCFDInt neles, const UCFDInt nvars,
                  const UCFDInt *restrict cell_ids,
                  const UCFDReal *restrict vol,
                  const UCFDReal *restrict du,
                  UCFDReal *restrict dup,
                  UCFDReal *restrict res)
{
    UCFDInt idx, ridx, kdx, offset;
    UCFDReal diff;

    OMPWrapper(ridx, kdx, diff)
    for (idx=0; idx<neles; ++idx)
    {
        ridx = cell_ids[idx];

        for (kdx=0; kdx<nvars; ++kdx) {
            offset = kdx + ridx*nvars;
            diff = du[offset] - dup[offset];
            res[idx + kdx*neles] = diff*diff*vol[idx];
            dup[offset] = du[offset];
        }
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_SubResidual(FlowSys sys, UCFDInt eidx)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;
    FlowElem *e  = &sys->eles[eidx];

    UCFDCall(rank_sub_residual(
        sys->nlocal, e->neles, sys->nvars, e->cell_ids,
        e->vol, blusgs->du, blusgs->dup, e->resid_out
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t
pre_blusgs(const UCFDInt nlocal, const UCFDInt nvars, const UCFDInt nfaces,
           const UCFDInt *restrict face_indptr, const UCFDInt *restrict face_slots,
           const UCFDInt8 *restrict face_sides, const UCFDReal *restrict face_area,
           const UCFDReal *restrict rcp_vol,
           UCFDReal *restrict diag, UCFDReal *restrict jmat)
{
    UCFDInt ridx, pos, row, col, slot, st, ed;
    UCFDInt8 side;
    const UCFDInt dim2 = nvars*nvars;
    UCFDReal dmat[dim2], fv, val;


    for (ridx=0; ridx<nlocal; ++ridx)
    {
        for (row=0; row<nvars; ++row) {
            for (col=0; col<nvars; ++col)
                dmat[col + row*nvars] = diag[col + row*nvars + ridx*dim2];
        }

        st = face_indptr[ridx];
        ed = face_indptr[ridx+1];
        for (pos=st; pos<ed; ++pos)
        {
            slot = face_slots[pos];
            side = face_sides[pos];
            fv = face_area[slot]*rcp_vol[ridx];

            for (row=0; row<nvars; ++row) {
                for (col=0; col<nvars; ++col) {
                    if (side == 1)
                        val = jmat[slot + col*nfaces + row*nfaces*nvars];
                    else
                        val = -jmat[slot + col*nfaces + row*nfaces*nvars + dim2*nfaces];
                    dmat[col + row*nvars] += val*fv;
                }
            }
        }
        ludcmp(nvars, dmat);

        for (row=0; row<nvars; ++row) {
            for (col=0; col<nvars; ++col)
                diag[col + row*nvars + ridx*dim2] = dmat[col + row*nvars];
        }
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_NSPrepare(FlowSys sys)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;

    UCFDCall(pre_blusgs(
        sys->nlocal, sys->nfvars, sys->nfaces,
        sys->rowptr, sys->slots, sys->sides, sys->face_area,
        sys->rcp_vol, blusgs->diag, blusgs->jmat
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_RANSPrepare(FlowSys sys)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;

    UCFDCall(pre_blusgs(
        sys->nlocal, sys->nturbvars, sys->nfaces,
        sys->rowptr, sys->slots, sys->sides, sys->face_area,
        sys->rcp_vol, blusgs->tdiag, blusgs->tjmat
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t
lower_sweep(const UCFDInt var0, const UCFDInt nv, const UCFDInt nlocal,
            const UCFDInt nvars, const UCFDInt ndims, const UCFDInt nfaces,
            const UCFDInt *restrict face_indptr, const UCFDInt *restrict face_neighbors,
            const UCFDInt8 *restrict face_sides, const UCFDInt *restrict face_slots,
            const UCFDReal *restrict face_area, const UCFDReal *restrict rcp_vol,
            const UCFDReal *restrict rhsb, const UCFDReal *restrict diag,
            const UCFDReal *restrict jmat, UCFDReal *restrict dub)
{
    UCFDInt ridx, row, col, pos, neib, slot;
    const UCFDInt dim2 = nv*nv;
    UCFDInt8 side;
    UCFDReal rhs[nv], dmat[dim2], fv, val, jval;

    for (ridx=0; ridx<nlocal; ++ridx)
    {
        for (row=0; row<nv; ++row) {
            rhs[row] = rhsb[var0+row + ridx*nvars];
            for (col=0; col<nv; ++col)
                dmat[col + nv*row] = diag[col + nv*row + ridx*dim2];
        }

        const UCFDInt pos_st = face_indptr[ridx];
        const UCFDInt pos_end = face_indptr[ridx+1];
        for (pos=pos_st; pos<pos_end; ++pos)
        {
            neib = face_neighbors[pos];
            slot = face_slots[pos];
            side = face_sides[pos];
            fv = face_area[slot]*rcp_vol[ridx];

            for (row=0; row<nv; ++row) {
                val = 0.0;
                for (col=0; col<nv; ++col) {
                    if (side == 1)
                        jval = jmat[slot + col*nfaces + row*nv*nfaces + dim2*nfaces];
                    else
                        jval = -jmat[slot + col*nfaces + row*nv*nfaces];
                    val += jval * dub[var0+col + neib*nvars];
                }
                rhs[row] -= val*fv;
            }
        }
        lusub(nv, dmat, rhs);
        for (row=0; row<nv; ++row)
            dub[var0+row + ridx*nvars] = rhs[row];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t
upper_sweep(const UCFDInt var0, const UCFDInt nv, const UCFDInt nlocal,
            const UCFDInt nvars, const UCFDInt ndims, const UCFDInt nfaces,
            const UCFDInt *restrict face_indptr, const UCFDInt *restrict face_neighbors,
            const UCFDInt8 *restrict face_sides, const UCFDInt *restrict face_slots,
            const UCFDReal *restrict face_area, const UCFDReal *restrict rcp_vol,
            const UCFDReal *restrict rhsb, UCFDReal *restrict diag,
            const UCFDReal *restrict jmat, UCFDReal *restrict dub)
{
    UCFDInt ridx, row, col, pos, neib, slot;
    const UCFDInt dim2 = nv*nv;
    UCFDInt8 side;
    UCFDReal rhs[nv], dmat[dim2], fv, val, jval;

    const UCFDInt i_begin = nlocal - 1;

    for (ridx=i_begin; ridx>-1; --ridx)
    {
        for (row=0; row<nv; ++row) {
            rhs[row] = rhsb[var0+row + ridx*nvars];
            for (col=0; col<nv; ++col)
                dmat[col + nv*row] = diag[col + nv*row + ridx*dim2];
        }

        const UCFDInt pos_st = face_indptr[ridx];
        const UCFDInt pos_end = face_indptr[ridx+1];
        for (pos=pos_st; pos<pos_end; ++pos)
        {
            neib = face_neighbors[pos];
            slot = face_slots[pos];
            side = face_sides[pos];
            fv = face_area[slot]*rcp_vol[ridx];

            for (row=0; row<nv; ++row) {
                val = 0.0;
                for (col=0; col<nv; ++col) {
                    if (side == 1)
                        jval = jmat[slot + col*nfaces + row*nv*nfaces + dim2*nfaces];
                    else
                        jval = -jmat[slot + col*nfaces + row*nv*nfaces];
                    val += jval * dub[var0+col + neib*nvars];
                }
                rhs[row] -= val*fv;
            }
        }
        lusub(nv, dmat, rhs);
        for (row=0; row<nv; ++row)
            dub[var0+row + ridx*nvars] = rhs[row];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_NSLowerSweep(FlowSys sys)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;

    UCFDCall(lower_sweep(
        0, sys->nfvars, sys->nlocal, sys->nvars, sys->ndims,
        sys->nfaces, sys->rowptr, sys->colidx, sys->sides,
        sys->slots, sys->face_area, sys->rcp_vol, blusgs->rhs,
        blusgs->diag, blusgs->jmat, blusgs->du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_RANSLowerSweep(FlowSys sys)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;
#if defined(DEBUG)
    UCFDCheckNull(blusgs->tjmat, "Turbulent jacobian matrix is not set\n");
#endif
    UCFDCall(lower_sweep(
        sys->nfvars, sys->nturbvars, sys->nlocal, sys->nvars, sys->ndims,
        sys->nfaces, sys->rowptr, sys->colidx, sys->sides,
        sys->slots, sys->face_area, sys->rcp_vol, blusgs->rhs,
        blusgs->tdiag, blusgs->tjmat, blusgs->du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_NSUpperSweep(FlowSys sys)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;

    UCFDCall(upper_sweep(
        0, sys->nfvars, sys->nlocal, sys->nvars, sys->ndims,
        sys->nfaces, sys->rowptr, sys->colidx, sys->sides,
        sys->slots, sys->face_area, sys->rcp_vol, blusgs->rhs,
        blusgs->diag, blusgs->jmat, blusgs->du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_RANSUpperSweep(FlowSys sys)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;
#if defined(DEBUG)
    UCFDCheckNull(blusgs->tjmat, "Turbulent jacobian matrix is not set\n");
#endif
    UCFDCall(upper_sweep(
        sys->nfvars, sys->nturbvars, sys->nlocal, sys->nvars, sys->ndims,
        sys->nfaces, sys->rowptr, sys->colidx, sys->sides,
        sys->slots, sys->face_area, sys->rcp_vol, blusgs->rhs,
        blusgs->tdiag, blusgs->tjmat, blusgs->du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_Reset(FlowSys sys)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;
    const UCFDInt nlocal = sys->nlocal, nvars = sys->nvars;

    memset(blusgs->du, 0, nlocal*nvars*sizeof(UCFDReal));
    memset(blusgs->dup, 0, nlocal*nvars*sizeof(UCFDReal));

    UCFDFunctionReturn(UCFD_SUCCESS);
}
