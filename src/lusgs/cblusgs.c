#include <string.h>
#include "flowsys.h"
#include "flux.h"
#include "smallmat.h"

/**
 * Overlapping function with BLU-SGS
 * - rank_blusgs_update
 * - rank_sub_residual
 * - rank_tblusgs_pack
 */

static ucfd_status_t
_rank_blusgs_diag_pack(const UCFDInt neles, const UCFDInt nlocal,
                       const UCFDInt nfvars, const UCFDReal a0,
                       const UCFDInt *restrict cell_ids,
                       const UCFDReal *restrict dt,
                       UCFDReal *restrict diag)
{
    UCFDInt idx, ridx, row, col;

    OMPWrapper(ridx, row, col)
    for (idx=0; idx<neles; ++idx) {
        ridx = cell_ids[idx];

        for (row=0; row<nfvars; ++row) {
            for (col=0; col<nfvars; ++col)
                diag[(size_t)(row*nfvars + col)*nlocal + ridx] = 0.0;
            diag[(size_t)(row*nfvars + row)*nlocal + ridx] = 1/dt[idx] + a0;
        }
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDCBLUSGS_Pack(UCFDSys sys, UCFDInt eidx, UCFDReal a0)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];

    UCFDCall(_rank_blusgs_diag_pack(
        e->neles, sys->nlocal, sys->nfvars, a0,
        e->cell_ids, e->dt, blusgs->flow_diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

/* turb pack - Use BLU-SGS function */
ucfd_status_t UCFDCBLUSGS_KWSST_Pack(UCFDSys sys, UCFDInt eidx, UCFDReal turb_factor, UCFDReal a0)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    srcjacobian dsrcf = kwsst_src_jacobian;

    UCFDCall(_rank_tblusgs_pack(
        e->neles, sys->nlocal, sys->nvars, sys->nfvars, a0, turb_factor,
        dsrcf, e->cell_ids, e->uptsb, e->dsrc, e->dt, blusgs->turb_diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDCBLUSGS_SA_Pack(UCFDSys sys, UCFDInt eidx, UCFDReal turb_factor, UCFDReal a0)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    srcjacobian dsrcf = sa_src_jacobian;

    UCFDCall(_rank_tblusgs_pack(
        e->neles, sys->nlocal, sys->nvars, sys->nfvars, a0, turb_factor,
        dsrcf, e->cell_ids, e->uptsb, e->dsrc, e->dt, blusgs->turb_diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

/* Update - Use BLU-SGS update function */
ucfd_status_t UCFDCBLUSGS_Update(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];

    UCFDCall(_rank_blusgs_update(
        e->neles, sys->nvars, sys->nlocal,
        e->cell_ids, blusgs->rank_du, e->uptsb
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDCBLUSGS_SubResidual(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];

    UCFDCall(_rank_sub_residual(
        e->neles, sys->nvars, sys->nlocal, e->cell_ids,
        e->vol, blusgs->rank_du, blusgs->rank_dup, e->resid_out
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t
_ele_pre_blusgs(const UCFDInt neles, const UCFDInt nv, const UCFDInt nface,
                const UCFDInt nlocal, const UCFDInt rank_nfaces,
                const UCFDInt *restrict cell_ids,
                const UCFDInt *restrict face_refs,
                const UCFDInt *restrict neighbors,
                const UCFDReal *restrict face_factors,
                UCFDReal *restrict diag,
                UCFDReal *restrict offdiag,
                const UCFDReal *restrict jmat)
{
    const UCFDInt dim2 = nv*nv;
    UCFDInt idx, ridx, fidx, row, col, face, slot;
    size_t ind0, ind1;
    UCFDReal dmat[dim2], fv, dval, oval;

    for (idx=0; idx<neles; ++idx)
    {
        ridx = cell_ids[idx];

        for (row=0; row<nv; ++row) {
            for (col=0; col<nv; ++col)
                dmat[col + row*nv] = diag[(size_t)(row*nv + col)*nlocal + ridx];
        }

        for (fidx=0; fidx<nface; ++fidx) {
            face = face_refs[idx + fidx*neles];
            UCFDBool has_neighbor = neighbors[idx + fidx*neles] >= 0;
            if (face > 0) slot = face - 1;
            else slot = -face - 1;
            fv = face_factors[idx + fidx*neles];

            for (row=0; row<nv; ++row) {
                for (col=0; col<nv; ++col) {
                    ind0 = slot + (size_t)col*rank_nfaces + (size_t)row*nv*rank_nfaces;
                    ind1 = ind0 + (size_t)dim2*rank_nfaces;
                    if (face > 0) dval = jmat[ind0];
                    else dval = -jmat[ind1];
                    dmat[col + row*nv] += dval*fv;

                    /* Colored sweeps use off-diagonal blocks */
                    if (has_neighbor) {
                        if (face > 0) oval = jmat[ind1];
                        else oval = -jmat[ind0];
                        offdiag[idx + fidx*neles + (size_t)col*nface*neles + (size_t)row*nv*nface*neles] = oval*fv;
                    }
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

ucfd_status_t UCFDCBLUSGS_NSPrepare(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];

    UCFDCall(_ele_pre_blusgs(
        e->neles, sys->nfvars, e->nface, sys->nlocal, sys->nfaces,
        e->cell_ids, e->face_refs, e->neighbors, e->face_factors,
        blusgs->flow_diag, e->flow_offdiag, blusgs->jmat
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDCBLUSGS_RANSPrepare(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];

    UCFDCall(_ele_pre_blusgs(
        e->neles, sys->nturbvars, e->nface, sys->nlocal, sys->nfaces,
        e->cell_ids, e->face_refs, e->neighbors, e->face_factors,
        blusgs->turb_diag, e->turb_offdiag, blusgs->tjmat
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t
_ele_colored_blusgs_sweep(const UCFDInt i_begin, const UCFDInt i_end,
                          const UCFDInt var0, const UCFDInt nv,
                          const UCFDInt neles, const UCFDInt nface,
                          const UCFDInt nlocal, const UCFDInt cell_offset,
                          const UCFDInt *restrict neighbors,
                          const UCFDReal *restrict diag ,
                          const UCFDReal *restrict offdiag,
                          const UCFDReal *restrict rhsb,
                          UCFDReal *restrict dub)
{
    const UCFDInt dim2 = nv*nv;
    UCFDInt idx, ridx, row, col, fidx, neib;
    UCFDReal rhs[nv], dmat[dim2], val;

    OMPWrapper(ridx, row, col, fidx, neib, rhs, dmat, val)
    for (idx=i_begin; idx<i_end; ++idx)
    {
        ridx = cell_offset + idx;

        for (row=0; row<nv; ++row) {
            rhs[row] = rhsb[idx + (size_t)(var0+row)*neles];
            for (col=0; col<nv; ++col)
                dmat[col + row*nv] = diag[(size_t)(row*nv + col)*nlocal + ridx];
        }

        for (fidx=0; fidx<nface; ++fidx) {
            neib = neighbors[idx + fidx*neles];
            if (neib >= 0) {
                for (row=0; row<nv; ++row) {
                    val = 0.0;
                    for (col=0; col<nv; ++col)
                        val += (
                            offdiag[idx + fidx*neles + (size_t)col*nface*neles + (size_t)row*nv*nface*neles]
                            * dub[neib + (var0+col)*nlocal]
                        );
                    rhs[row] -= val;
                }
            }
        }

        block_matvec(nv, dmat, rhs);

        for (row=0; row<nv; ++row)
            dub[ridx + (var0+row)*nlocal] = rhs[row];
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDCBLUSGS_NSLowerSweep(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    const UCFDInt ncolors = e->ncolors;

    for (UCFDInt i=0; i<ncolors; ++i)
        UCFDCall(_ele_colored_blusgs_sweep(
            e->color_offset[i], e->color_offset[i+1], 0, sys->nfvars,
            e->neles, e->nface, sys->nlocal, e->cell_offset, e->neighbors,
            blusgs->flow_diag, e->flow_offdiag, e->rhs, blusgs->rank_du
        ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDCBLUSGS_NSUpperSweep(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    const UCFDInt ncolors = e->ncolors;

    for (UCFDInt i=ncolors-2; i>=0; --i)
        UCFDCall(_ele_colored_blusgs_sweep(
            e->color_offset[i], e->color_offset[i+1], 0, sys->nfvars,
            e->neles, e->nface, sys->nlocal, e->cell_offset, e->neighbors,
            blusgs->flow_diag, e->flow_offdiag, e->rhs, blusgs->rank_du
        ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDCBLUSGS_RANSLowerSweep(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    const UCFDInt ncolors = e->ncolors;

    for (UCFDInt i=0; i<ncolors; ++i)
        UCFDCall(_ele_colored_blusgs_sweep(
            e->color_offset[i], e->color_offset[i+1], sys->nfvars, sys->nturbvars,
            e->neles, e->nface, sys->nlocal, e->cell_offset, e->neighbors,
            blusgs->turb_diag, e->turb_offdiag, e->rhs, blusgs->rank_du
        ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDCBLUSGS_RANSUpperSweep(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    const UCFDInt ncolors = e->ncolors;

    for (UCFDInt i=ncolors-2; i>=0; --i)
        UCFDCall(_ele_colored_blusgs_sweep(
            e->color_offset[i], e->color_offset[i+1], sys->nfvars, sys->nturbvars,
            e->neles, e->nface, sys->nlocal, e->cell_offset, e->neighbors,
            blusgs->turb_diag, e->turb_offdiag, e->rhs, blusgs->rank_du
        ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDCBLUSGS_Reset(UCFDSys sys)
{
    return UCFDBLUSGS_Reset(sys);
}


ucfd_status_t UCFDColoredSysSetBLUSGS(UCFDSys *sys,
                                      UCFDReal *rank_du, UCFDReal *rank_dup,
                                      UCFDReal *flow_diag, UCFDReal *jmat,
                                      UCFDReal *turb_diag, UCFDReal *tjmat)
{
    UCFDSys s           = *sys;
    PFlowSys *fsys      = (PFlowSys *)s->flowsys;

    fsys->rank_du        = rank_du;
    fsys->rank_dup       = rank_dup;
    fsys->flow_diag      = flow_diag;
    fsys->jmat           = jmat;
    fsys->turb_diag      = turb_diag;
    fsys->tjmat          = tjmat;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDColoredEleSetBLUSGS(UCFDSys *sys, UCFDInt eidx,
                                      UCFDInt neles, UCFDInt nface,
                                      UCFDInt ncolors, UCFDInt cell_offset,
                                      UCFDInt *nblocks, UCFDInt *color_offset,
                                      UCFDInt *cell_ids, UCFDInt *face_refs,
                                      UCFDInt *neighbors, UCFDReal *face_factors,
                                      UCFDReal *uptsb, UCFDReal *rhs,
                                      UCFDReal *dt, UCFDReal *dsrc,
                                      UCFDReal *vol, UCFDReal *resid_out,
                                      UCFDReal *flow_offdiag, UCFDReal *turb_offdiag)
{
    UCFDSys s           = *sys;
    PFlowEle *eles      = (PFlowEle *)s->eles;
    PFlowEle *e         = &eles[eidx];

    e->neles            = neles;
    e->nface            = nface;
    e->cell_ids         = cell_ids;
    e->uptsb            = uptsb;
    e->rhs              = rhs;
    e->dt               = dt;
    e->dsrc             = dsrc;
    e->vol              = vol;
    e->resid_out        = resid_out;
    e->face_refs        = face_refs;
    e->neighbors        = neighbors;
    e->face_factors     = face_factors;
    e->flow_offdiag     = flow_offdiag;
    e->turb_offdiag     = turb_offdiag;

    /* Coloring properties */
    e->nblocks          = nblocks;
    e->ncolors          = ncolors;
    e->cell_offset      = cell_offset;
    e->color_offset     = color_offset;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
