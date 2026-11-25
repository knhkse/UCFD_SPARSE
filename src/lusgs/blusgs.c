#include <string.h>
#include "flowsys.h"
#include "flux.h"
#include "smallmat.h"


/**
 * ----------------- array structure -------------
 * rhs, du, dup : [nvars, nlocal]
 * flow_diag : [nlocal, nfvars, nfvars]
 * turb_diag : [nlocal, nturbvars, nturbvars]
 */

/**
 * Pack & update kernels => per-element execution
 */
static ucfd_status_t
_rank_blusgs_pack(const UCFDInt neles, const UCFDInt nlocal,
                  const UCFDInt nvars, const UCFDInt nfvars, 
                  const UCFDReal a0,
                  const UCFDInt *restrict cell_ids,
                  const UCFDReal *restrict rhs,
                  const UCFDReal *restrict dt,
                  UCFDReal *restrict rank_rhs,
                  UCFDReal *restrict flow_diag)
{
    UCFDInt idx, ridx, kdx, row, col;

    for (idx=0; idx<neles; ++idx) {
        ridx = cell_ids[idx];
        for (kdx=0; kdx<nvars; ++kdx)
            // rank_rhs[ridx, kdx] = rhs[kdx, idx]
            rank_rhs[ridx + kdx*nlocal] = rhs[idx + kdx*neles];
        
        for (row=0; row<nfvars; ++row) {
            for (col=0; col<nfvars; ++col)
                flow_diag[(size_t)(row*nfvars + col)*nlocal + ridx] = 0.0;
            flow_diag[(size_t)(row*nfvars + row)*nlocal + ridx] = 1/dt[idx] + a0;
        }
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_Pack(UCFDSys sys, UCFDInt eidx, UCFDReal a0)
{
    FlowSys *blusgs = (FlowSys *)sys->flowsys;
    FlowEle *e  = &((FlowEle *)sys->eles)[eidx];
    
    UCFDCall(_rank_blusgs_pack(
        e->neles, sys->nlocal, sys->nvars, sys->nfvars, a0,
        e->cell_ids, e->rhs, e->dt, blusgs->rank_rhs, blusgs->flow_diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


ucfd_status_t
_rank_tblusgs_pack(const UCFDInt neles, const UCFDInt nlocal,
                   const UCFDInt nvars, const UCFDInt nfvars,
                   const UCFDReal a0, const UCFDReal factor,
                   srcjacobian dsrcf,
                   const UCFDInt *restrict cell_ids,
                   const UCFDReal *restrict uptsb,
                   const UCFDReal *restrict dsrc,
                   const UCFDReal *restrict dt,
                   UCFDReal *restrict turb_diag)
{
    UCFDInt idx, ridx, kdx, row, col;
    const UCFDInt nturbvars = nvars - nfvars;
    UCFDReal u[nvars], ds[nvars], dmat[nturbvars*nturbvars];

    for (idx=0; idx<neles; ++idx)
    {
        ridx = cell_ids[idx];
        for (row=0; row<nturbvars; ++row) {
            for (col=0; col<nturbvars; ++col)
                dmat[col+row*nturbvars] = 0.0;
        }

        // Prepare cell uf and dsrc
        for (kdx=0; kdx<nvars; ++kdx) {
            u[kdx] = uptsb[idx + kdx*neles];
            ds[kdx] = dsrc[idx + kdx*neles];
        }
        dsrcf(nvars, nturbvars, u, dmat, ds);

        for (row=0; row<nturbvars; ++row) {
            for (col=0; col<nturbvars; ++col)
                turb_diag[(size_t)(row*nturbvars + col)*nlocal + ridx] = dmat[col+row*nturbvars];
            turb_diag[(size_t)(row*nturbvars + row)*nlocal + ridx] += 1.0/(dt[idx]*factor) + a0;
        }
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_KWSST_Pack(UCFDSys sys, UCFDInt eidx, UCFDReal turb_factor, UCFDReal a0)
{
    FlowSys *blusgs = (FlowSys *)sys->flowsys;
    FlowEle *e  = &((FlowEle *)sys->eles)[eidx];
    srcjacobian dsrcf = kwsst_src_jacobian;

    UCFDCall(_rank_tblusgs_pack(
        e->neles, sys->nlocal, sys->nvars, sys->nfvars, a0, turb_factor,
        dsrcf, e->cell_ids, e->uptsb, e->dsrc, e->dt, blusgs->turb_diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_SA_Pack(UCFDSys sys, UCFDInt eidx, UCFDReal turb_factor, UCFDReal a0)
{
    FlowSys *blusgs = (FlowSys *)sys->flowsys;
    FlowEle *e  = &((FlowEle *)sys->eles)[eidx];
    srcjacobian dsrcf = sa_src_jacobian;

    UCFDCall(_rank_tblusgs_pack(
        e->neles, sys->nlocal, sys->nvars, sys->nfvars, a0, turb_factor,
        dsrcf, e->cell_ids, e->uptsb, e->dsrc, e->dt, blusgs->turb_diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


/**
 * Update
 */
ucfd_status_t
_rank_blusgs_update(const UCFDInt neles,
                    const UCFDInt nvars, const UCFDInt nlocal,
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
            upts[idx + kdx*neles] += rank_du[ridx + kdx*nlocal];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_Update(UCFDSys sys, UCFDInt eidx)
{
    FlowSys *blusgs = (FlowSys *)sys->flowsys;
    FlowEle *e = &((FlowEle *)sys->eles)[eidx];

    UCFDCall(_rank_blusgs_update(
        e->neles, sys->nvars, sys->nlocal,
        e->cell_ids, blusgs->rank_du, e->uptsb
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


/**
 * Inner iteration residual
 */
ucfd_status_t
_rank_sub_residual(const UCFDInt neles,
                   const UCFDInt nvars, const UCFDInt nlocal,
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
            offset = ridx + kdx*nlocal;
            diff = du[offset] - dup[offset];
            res[idx + kdx*neles] = diff*diff*vol[idx];
            dup[offset] = du[offset];
        }
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_SubResidual(UCFDSys sys, UCFDInt eidx)
{
    FlowSys *blusgs = (FlowSys *)sys->flowsys;
    FlowEle *e  = &((FlowEle *)sys->eles)[eidx];

    UCFDCall(_rank_sub_residual(
        e->neles, sys->nvars, sys->nlocal, e->cell_ids,
        e->vol, blusgs->rank_du, blusgs->rank_dup, e->resid_out
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t
_pre_blusgs(const UCFDInt nlocal, const UCFDInt nv, const UCFDInt nfaces,
            const UCFDInt *restrict face_indptr, const UCFDInt *restrict face_slots,
            const UCFDInt8 *restrict face_sides, const UCFDReal *restrict face_area,
            const UCFDReal *restrict rcp_vol, UCFDReal *restrict diag, UCFDReal *restrict jmat)
{
    UCFDInt ridx, pos, row, col, slot, st, ed;
    UCFDInt8 side;
    const UCFDInt dim2 = nv*nv;
    UCFDReal dmat[dim2], fv, val;

    for (ridx=0; ridx<nlocal; ++ridx)
    {
        for (row=0; row<nv; ++row) {
            for (col=0; col<nv; ++col)
                dmat[col + row*nv] = diag[(size_t)(row*nv + col)*nlocal + ridx];
        }

        st = face_indptr[ridx];
        ed = face_indptr[ridx+1];
        for (pos=st; pos<ed; ++pos)
        {
            slot = face_slots[pos];
            side = face_sides[pos];
            fv = face_area[slot]*rcp_vol[ridx];

            for (row=0; row<nv; ++row) {
                for (col=0; col<nv; ++col) {
                    if (side == 1)
                        val = jmat[slot + col*nfaces + row*nfaces*nv];
                    else
                        val = -jmat[slot + col*nfaces + row*nfaces*nv + dim2*nfaces];
                    dmat[col + row*nv] += val*fv;
                }
            }
        }
        pivoted_inverse(nv, dmat, PIVOT_TOL);

        for (row=0; row<nv; ++row) {
            for (col=0; col<nv; ++col)
                diag[(size_t)(row*nv + col)*nlocal + ridx] = dmat[col + row*nv];
        }
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_NSPrepare(UCFDSys sys)
{
    FlowSys *blusgs = (FlowSys *)sys->flowsys;

    UCFDCall(_pre_blusgs(
        sys->nlocal, sys->nfvars, sys->nfaces,
        blusgs->face_indptr, blusgs->face_slots, blusgs->face_sides,
        blusgs->face_area, blusgs->rcp_vol, blusgs->flow_diag, blusgs->jmat
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_RANSPrepare(UCFDSys sys)
{
    FlowSys *blusgs = (FlowSys *)sys->flowsys;

    UCFDCall(_pre_blusgs(
        sys->nlocal, sys->nturbvars, sys->nfaces,
        blusgs->face_indptr, blusgs->face_slots, blusgs->face_sides,
        blusgs->face_area, blusgs->rcp_vol, blusgs->turb_diag, blusgs->tjmat
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t
_lower_sweep(const UCFDInt var0, const UCFDInt nv,
             const UCFDInt nlocal, const UCFDInt nfaces,
             const UCFDInt *restrict face_indptr, const UCFDInt *restrict face_neighbors,
             const UCFDInt8 *restrict face_sides, const UCFDInt *restrict face_slots,
             const UCFDReal *restrict face_area, const UCFDReal *restrict rcp_vol,
             const UCFDReal *restrict rank_rhs, const UCFDReal *restrict diag,
             const UCFDReal *restrict jmat, UCFDReal *restrict rank_du)
{
    UCFDInt ridx, row, col, pos, neib, slot;
    const UCFDInt dim2 = nv*nv;
    UCFDInt8 side;
    UCFDReal rhs[nv], dmat[dim2], fv, val, jval;

    for (ridx=0; ridx<nlocal; ++ridx)
    {
        for (row=0; row<nv; ++row) {
            rhs[row] = rank_rhs[ridx + (var0+row)*nlocal];
            for (col=0; col<nv; ++col)
                dmat[col + nv*row] = diag[(size_t)(row*nv + col)*nlocal + ridx];
        }

        const UCFDInt pos_st = face_indptr[ridx];
        const UCFDInt pos_end = face_indptr[ridx+1];
        for (pos=pos_st; pos<pos_end; ++pos)
        {
            neib = face_neighbors[pos];
            if (neib == ridx) continue;

            slot = face_slots[pos];
            side = face_sides[pos];
            fv = face_area[slot]*rcp_vol[ridx];

            for (row=0; row<nv; ++row) {
                val = 0.0;
                for (col=0; col<nv; ++col) {
                    if (side == 1)
                        jval = jmat[(size_t)(dim2+row*nv+col)*nfaces + slot];
                    else
                        jval = -jmat[(size_t)(row*nv+col)*nfaces + slot];
                    val += jval * rank_du[neib + (var0+col)*nlocal];
                }
                rhs[row] -= val*fv;
            }
        }
        block_matvec(nv, dmat, rhs);
        for (row=0; row<nv; ++row)
            rank_du[ridx + (var0+row)*nlocal] = rhs[row];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t
_upper_sweep(const UCFDInt var0, const UCFDInt nv,
             const UCFDInt nlocal, const UCFDInt nfaces,
             const UCFDInt *restrict face_indptr, const UCFDInt *restrict face_neighbors,
             const UCFDInt8 *restrict face_sides, const UCFDInt *restrict face_slots,
             const UCFDReal *restrict face_area, const UCFDReal *restrict rcp_vol,
             const UCFDReal *restrict rank_rhs, UCFDReal *restrict diag,
             const UCFDReal *restrict jmat, UCFDReal *restrict rank_du)
{
    UCFDInt ridx, row, col, pos, neib, slot;
    const UCFDInt dim2 = nv*nv;
    UCFDInt8 side;
    UCFDReal rhs[nv], dmat[dim2], fv, val, jval;

    const UCFDInt i_begin = nlocal - 1;

    for (ridx=i_begin; ridx>-1; --ridx)
    {
        for (row=0; row<nv; ++row) {
            rhs[row] = rank_rhs[ridx + (var0+row)*nlocal];
            for (col=0; col<nv; ++col)
                dmat[col + nv*row] = diag[(size_t)(row*nv+col)*nlocal + ridx];
        }

        const UCFDInt pos_st = face_indptr[ridx];
        const UCFDInt pos_end = face_indptr[ridx+1];
        for (pos=pos_st; pos<pos_end; ++pos)
        {
            neib = face_neighbors[pos];
            if (neib == ridx) continue;

            slot = face_slots[pos];
            side = face_sides[pos];
            fv = face_area[slot]*rcp_vol[ridx];

            for (row=0; row<nv; ++row) {
                val = 0.0;
                for (col=0; col<nv; ++col) {
                    if (side == 1)
                        jval = jmat[(size_t)(dim2+row*nv+col)*nfaces + slot];
                    else
                        jval = -jmat[(size_t)(row*nv+col)*nfaces + slot];
                    val += jval * rank_du[neib + (var0+col)*nlocal];
                }
                rhs[row] -= val*fv;
            }
        }
        block_matvec(nv, dmat, rhs);
        for (row=0; row<nv; ++row)
            rank_du[ridx + (var0+row)*nlocal] = rhs[row];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_NSLowerSweep(UCFDSys sys)
{
    FlowSys *blusgs = (FlowSys *)sys->flowsys;

    UCFDCall(_lower_sweep(
        0, sys->nfvars, sys->nlocal, sys->nfaces,
        blusgs->face_indptr, blusgs->face_neighbors, blusgs->face_sides,
        blusgs->face_slots, blusgs->face_area, blusgs->rcp_vol, blusgs->rank_rhs,
        blusgs->flow_diag, blusgs->jmat, blusgs->rank_du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_RANSLowerSweep(UCFDSys sys)
{
    FlowSys *blusgs = (FlowSys *)sys->flowsys;
#if defined(DEBUG)
    UCFDCheckNull(blusgs->tjmat, "Turbulent jacobian matrix is not set\n");
#endif
    UCFDCall(_lower_sweep(
        sys->nfvars, sys->nturbvars, sys->nlocal, sys->nfaces,
        blusgs->face_indptr, blusgs->face_neighbors, blusgs->face_sides,
        blusgs->face_slots, blusgs->face_area, blusgs->rcp_vol, blusgs->rank_rhs,
        blusgs->turb_diag, blusgs->tjmat, blusgs->rank_du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_NSUpperSweep(UCFDSys sys)
{
    FlowSys *blusgs = (FlowSys *)sys->flowsys;

    UCFDCall(_upper_sweep(
        0, sys->nfvars, sys->nlocal, sys->nfaces,
        blusgs->face_indptr, blusgs->face_neighbors, blusgs->face_sides,
        blusgs->face_slots, blusgs->face_area, blusgs->rcp_vol, blusgs->rank_rhs,
        blusgs->flow_diag, blusgs->jmat, blusgs->rank_du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_RANSUpperSweep(UCFDSys sys)
{
    FlowSys *blusgs = (FlowSys *)sys->flowsys;
#if defined(DEBUG)
    UCFDCheckNull(blusgs->tjmat, "Turbulent jacobian matrix is not set\n");
#endif
    UCFDCall(_upper_sweep(
        sys->nfvars, sys->nturbvars, sys->nlocal, sys->nfaces,
        blusgs->face_indptr, blusgs->face_neighbors, blusgs->face_sides,
        blusgs->face_slots, blusgs->face_area, blusgs->rcp_vol, blusgs->rank_rhs,
        blusgs->turb_diag, blusgs->tjmat, blusgs->rank_du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDBLUSGS_Reset(UCFDSys sys)
{
    FlowSys *blusgs = (FlowSys *)sys->flowsys;
    const UCFDInt nlocal = sys->nlocal, nvars = sys->nvars;

    memset(blusgs->rank_du, 0, nlocal*nvars*sizeof(UCFDReal));
    memset(blusgs->rank_dup, 0, nlocal*nvars*sizeof(UCFDReal));

    UCFDFunctionReturn(UCFD_SUCCESS);
}


ucfd_status_t UCFDSysSetBLUSGS(UCFDSys *sys,
                               UCFDReal *rank_du, UCFDReal *rank_dup, UCFDReal *rank_rhs,
                               UCFDReal *flow_diag, UCFDReal *jmat,
                               UCFDReal *turb_diag, UCFDReal *tjmat)
{
    UCFDSys s           = *sys;
    FlowSys *fsys       = (FlowSys *)s->flowsys;

    fsys->rank_du       = rank_du;
    fsys->rank_dup      = rank_dup;
    fsys->rank_rhs      = rank_rhs;
    fsys->flow_diag     = flow_diag;
    fsys->jmat          = jmat;
    fsys->turb_diag     = turb_diag;
    fsys->tjmat         = tjmat;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDEleSetBLUSGS(UCFDSys *sys, UCFDInt eidx,
                              UCFDInt neles, UCFDInt nface, UCFDInt *cell_ids,
                              UCFDReal *uptsb, UCFDReal *rhs,
                              UCFDReal *dt, UCFDReal *dsrc,
                              UCFDReal *vol, UCFDReal *resid_out)
{
    UCFDSys s           = *sys;
    FlowEle *eles       = (FlowEle *)s->eles;
    FlowEle *e          = &eles[eidx];

    e->neles            = neles;
    e->nface            = nface;
    e->cell_ids         = cell_ids;
    e->uptsb            = uptsb;
    e->rhs              = rhs;
    e->dt               = dt;
    e->dsrc             = dsrc;
    e->vol              = vol;
    e->resid_out        = resid_out;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
