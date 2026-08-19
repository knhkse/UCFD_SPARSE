#include "flowsys.h"
#include "flux.h"


/**
 * Pack & update kernels => per-element execution
 */
static ucfd_status_t
lusgs_pack(const UCFDInt nlocal, const UCFDInt neles, const UCFDInt nvars,
           const UCFDInt nfvars, const UCFDReal turb_factor, const UCFDReal a0,
           const UCFDInt *restrict cell_ids,
           const UCFDReal *restrict upts,
           const UCFDReal *restrict rhs,
           const UCFDReal *restrict dt,
           const UCFDReal *restrict dsrc,
           UCFDReal *restrict rank_u,
           UCFDReal *restrict rank_du,
           UCFDReal *restrict rank_diag)
{
    UCFDInt idx, ridx, kdx, rank_idx, ele_idx;
    UCFDReal factor;

    OMPWrapper(ridx, kdx, rank_idx, ele_idx, factor)
    for (idx=0; idx<neles; ++idx)
    {
        ridx = cell_ids[idx];
        for (kdx=0; kdx<nvars; ++kdx) {
            rank_idx = ridx + kdx*nlocal;
            ele_idx = idx + kdx*neles;
            factor = kdx < nfvars ? 1.0 : turb_factor;
            rank_u[rank_idx] = upts[ele_idx];
            rank_du[rank_idx] = rhs[ele_idx];
            rank_diag[rank_idx] = 1.0/(dt[idx]*factor) + a0 + dsrc[ele_idx];
        }
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_Pack(FlowSys sys, UCFDInt eidx,
                             UCFDReal a0, UCFDReal turb_factor)
{
    LUSGSSys *lusgs = (LUSGSSys *)sys->data;
    FlowElem *e  = &(sys->eles[eidx]);

    UCFDCall(lusgs_pack(
        sys->nlocal, e->neles, sys->nvars, sys->nfvars, turb_factor, a0,
        e->cell_ids, e->uptsb, e->rhs, e->dt, e->dsrc,
        sys->u, sys->du, lusgs->diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t
rank_lusgs_update(const UCFDInt nlocal, const UCFDInt neles, const UCFDInt nvars,
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

ucfd_status_t UCFDLUSGS_Update(FlowSys sys, UCFDInt eidx)
{
    FlowElem *e  = &sys->eles[eidx];

    UCFDCall(rank_lusgs_update(
        sys->nlocal, e->neles, sys->nvars,
        e->cell_ids, sys->du, e->uptsb
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

/**
 * Preparation => Compute diagonal element
 */
static ucfd_status_t
pre_lusgs(const UCFDInt nlocal, const UCFDInt nv0, const UCFDInt nv1,
          const UCFDReal kappa,
          const UCFDInt *restrict face_indptr, const UCFDInt *restrict face_slots,
          const UCFDReal *restrict face_area, const UCFDReal *restrict rcp_vol,
          const UCFDReal *restrict lambdaf, UCFDReal *restrict diag)
{
    UCFDInt ridx, pos, slot, kdx, st, ed;
    UCFDReal lamf, spectral_diag;

    OMPWrapper(st, ed, spectral_diag, pos, slot, lamf, kdx)
    for (ridx=0; ridx<nlocal; ++ridx)
    {
        st = face_indptr[ridx];
        ed = face_indptr[ridx+1];
        spectral_diag = 0.0;
        for (pos=st; pos<ed; ++pos)
        {
            slot = face_slots[pos];
            lamf = lambdaf[slot]*kappa;
            spectral_diag += 0.5*lamf*face_area[slot]*rcp_vol[ridx];
        }
        for (kdx=nv0; kdx<nv1; ++kdx)
            diag[ridx+kdx*nlocal] += spectral_diag;
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_NSPrepare(FlowSys sys, UCFDReal kappa)
{
    LUSGSSys *lusgs = (LUSGSSys *)sys->data;

    UCFDCall(pre_lusgs(
        sys->nlocal, 0, sys->nfvars, kappa,
        sys->rowptr, sys->slots, sys->face_area, sys->rcp_vol,
        lusgs->fspr, lusgs->diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_RANSPrepare(FlowSys sys, UCFDReal kappa)
{
    LUSGSSys *lusgs = (LUSGSSys *)sys->data;

    UCFDCall(pre_lusgs(
        sys->nlocal, sys->nfvars, sys->nvars, kappa,
        sys->rowptr, sys->slots, sys->face_area, sys->rcp_vol,
        lusgs->tfspr, lusgs->diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


/**
 * Sweep functions
 */
static inline void diff_flux(UCFDInt nvars, UCFDInt nfvars, UCFDInt nturbvars,
                             UCFDInt ndims, UCFDInt dnv,
                             fluxfunc fluxf,
                             UCFDReal *restrict u, UCFDReal *restrict du, 
                             UCFDReal *restrict df, UCFDReal *restrict nf)
{
    UCFDInt i;
    UCFDReal f[dnv];
    for (i=0; i<nvars; ++i) du[i] += u[i];
    fluxf(nfvars, nturbvars, ndims, u, nf, f);
    fluxf(nfvars, nturbvars, ndims, du, nf, df);
    for (i=0; i<dnv; ++i) df[i] -= f[i];
}

static ucfd_status_t
lower_sweep(const UCFDInt nv0, const UCFDInt nv1, const UCFDReal kappa,
            fluxfunc fluxf, const UCFDReal *restrict lambdaf,
            const UCFDInt nlocal, const UCFDInt nvars, const UCFDInt nfvars,
            const UCFDInt ndims, const UCFDInt nfaces,
            const UCFDInt *restrict face_indptr, const UCFDInt *restrict face_neighbors,
            const UCFDInt8 *restrict face_sides, const UCFDInt *restrict face_slots,
            const UCFDReal *restrict face_area, const UCFDReal *restrict face_normal,
            const UCFDReal *restrict rcp_vol, const UCFDReal *restrict diag,
            const UCFDReal *restrict rank_u, UCFDReal *restrict rank_du)
{
    UCFDInt ridx, kdx, pos, neib, slot;
    UCFDInt8 side;
    const UCFDInt dnv = nv1 - nv0, ntvars = nvars - nfvars;
    UCFDReal du[nvars], dfj[dnv], df[dnv];
    UCFDReal u[nvars], nf[ndims];
    UCFDReal fv;

    for (ridx=0; ridx<nlocal; ++ridx)
    {
        for (kdx=0; kdx<dnv; ++kdx) df[kdx] = 0.0;

        const UCFDInt pos_st = face_indptr[ridx];
        const UCFDInt pos_end = face_indptr[ridx+1];
        for (pos=pos_st; pos<pos_end; ++pos)
        {
            neib = face_neighbors[pos];
            if (neib < ridx) {
                slot = face_slots[pos];
                for (kdx=0; kdx<nvars; ++kdx) u[kdx] = rank_u[neib+kdx*nlocal];
                side = face_sides[pos];
                for (kdx=0; kdx<ndims; ++kdx) nf[kdx] = side*face_normal[slot+kdx*nfaces];
                fv = face_area[slot] * rcp_vol[ridx];

                for (kdx=0; kdx<nvars; ++kdx) du[kdx] = 0.0;
                for (kdx=nv0; kdx<nv1; ++kdx) du[kdx] = rank_du[neib + kdx*nlocal];
                diff_flux(nvars, nfvars, ntvars, ndims, dnv, fluxf, u, du, dfj, nf);

                for (kdx=0; kdx<dnv; ++kdx)
                    df[kdx] += (dfj[kdx] - kappa*lambdaf[slot]*rank_du[neib+(kdx+nv0)*nlocal])*fv;
            }
        }
        for (kdx=0; kdx<dnv; ++kdx) {
            const UCFDInt ldx = ridx + (kdx+nv0)*nlocal;
            rank_du[ldx] = (rank_du[ldx] - 0.5*df[kdx])/diag[ldx];
        }
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t
upper_sweep(const UCFDInt nv0, const UCFDInt nv1, const UCFDReal kappa,
            fluxfunc fluxf, const UCFDReal *restrict lambdaf,
            const UCFDInt nlocal, const UCFDInt nvars, const UCFDInt nfvars,
            const UCFDInt ndims, const UCFDInt nfaces,
            const UCFDInt *restrict face_indptr, const UCFDInt *restrict face_neighbors,
            const UCFDInt8 *restrict face_sides, const UCFDInt *restrict face_slots,
            const UCFDReal *restrict face_area, const UCFDReal *restrict face_normal,
            const UCFDReal *restrict rcp_vol, const UCFDReal *restrict diag,
            const UCFDReal *restrict rank_u, UCFDReal *restrict rank_du)
{
    UCFDInt ridx, kdx, pos, neib, slot;
    UCFDInt8 side;
    const UCFDInt dnv = nv1 - nv0, ntvars = nvars - nfvars;
    UCFDReal du[nvars], dfj[dnv], df[dnv];
    UCFDReal u[nvars], nf[ndims];
    UCFDReal fv;

    const UCFDInt i_begin = nlocal-1;

    for (ridx=i_begin; ridx>-1; --ridx)
    {
        for (kdx=0; kdx<dnv; ++kdx) df[kdx] = 0.0;

        const UCFDInt pos_st = face_indptr[ridx];
        const UCFDInt pos_end = face_indptr[ridx+1];
        for (pos=pos_st; pos<pos_end; ++pos)
        {
            neib = face_neighbors[pos];
            if (neib > ridx) {
                slot = face_slots[pos];
                for (kdx=0; kdx<nvars; ++kdx) u[kdx] = rank_u[neib+kdx*nlocal];
                side = face_sides[pos];
                for (kdx=0; kdx<ndims; ++kdx) nf[kdx] = side*face_normal[slot+kdx*nfaces];
                fv = face_area[slot] * rcp_vol[ridx];

                for (kdx=0; kdx<nvars; ++kdx) du[kdx] = 0.0;
                for (kdx=nv0; kdx<nv1; ++kdx) du[kdx] = rank_du[neib + kdx*nlocal];
                diff_flux(nvars, nfvars, ntvars, ndims, dnv, fluxf, u, du, dfj, nf);

                for (kdx=0; kdx<dnv; ++kdx)
                    df[kdx] += (dfj[kdx] - kappa*lambdaf[slot]*rank_du[neib+(kdx+nv0)*nlocal])*fv;
            }
        }
        for (kdx=0; kdx<dnv; ++kdx) {
            const UCFDInt ldx = ridx + (kdx+nv0)*nlocal;
            rank_du[ldx] -= (0.5*df[kdx] / diag[ldx]);
        }
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_NSLowerSweep(FlowSys sys, UCFDReal kappa)
{
    LUSGSSys *lusgs = (LUSGSSys *)sys->data;
    fluxfunc f = ns_flux_container;
    
    UCFDCall(lower_sweep(
        0, sys->nfvars, kappa, f, lusgs->fspr,
        sys->nlocal, sys->nvars, sys->nfvars, sys->ndims, sys->nfaces,
        sys->rowptr, sys->colidx, sys->sides, sys->slots, sys->face_area,
        sys->face_normal, sys->rcp_vol, lusgs->diag,
        sys->u, sys->du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_RANSLowerSweep(FlowSys sys, UCFDReal kappa)
{
#if defined(DEBUG)
    UCFDCheckNull(sys->tfspr, "Turbulent spectral radius is not set\n");
#endif
    LUSGSSys *lusgs = (LUSGSSys *)sys->data;
    fluxfunc f = rans_flux_container;
    
    UCFDCall(lower_sweep(
        sys->nfvars, sys->nvars, kappa, f, lusgs->tfspr,
        sys->nlocal, sys->nvars, sys->nfvars, sys->ndims, sys->nfaces,
        sys->rowptr, sys->colidx, sys->sides, sys->slots, sys->face_area,
        sys->face_normal, sys->rcp_vol, lusgs->diag,
        sys->u, sys->du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_NSUpperSweep(FlowSys sys, UCFDReal kappa)
{
    LUSGSSys *lusgs = (LUSGSSys *)sys->data;
    fluxfunc f = ns_flux_container;
    
    UCFDCall(upper_sweep(
        0, sys->nfvars, kappa, f, lusgs->fspr,
        sys->nlocal, sys->nvars, sys->nfvars, sys->ndims, sys->nfaces,
        sys->rowptr, sys->colidx, sys->sides, sys->slots, sys->face_area,
        sys->face_normal, sys->rcp_vol, lusgs->diag,
        sys->u, sys->du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_RANSUpperSweep(FlowSys sys, UCFDReal kappa)
{
#if defined(DEBUG)
    UCFDCheckNull(sys->tfspr, "Turbulent spectral radius is not set\n");
#endif
    LUSGSSys *lusgs = (LUSGSSys *)sys->data;
    fluxfunc f = rans_flux_container;
    
    UCFDCall(upper_sweep(
        sys->nfvars, sys->nvars, kappa, f, lusgs->tfspr,
        sys->nlocal, sys->nvars, sys->nfvars, sys->ndims, sys->nfaces,
        sys->rowptr, sys->colidx, sys->sides, sys->slots, sys->face_area,
        sys->face_normal, sys->rcp_vol, lusgs->diag,
        sys->u, sys->du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}
