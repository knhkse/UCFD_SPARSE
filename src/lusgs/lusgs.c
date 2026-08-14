#include "flowsys.h"
#include "flux.h"


// ! Test functions
void TestPack(FlowSys sys, UCFDInt eidx)
{
    FlowElem *e  = &sys->eles[eidx];
    const UCFDInt nvars = sys->nvars;
    const UCFDInt *cell_ids = e->cell_ids;
    UCFDReal *upts = e->uptsb, *rhs = e->rhs;
    UCFDReal *rank_u = sys->u, *rank_du = sys->du;
    const UCFDInt neles = e->neles, nlocal = sys->nlocal;

    for (int idx=0; idx<neles; ++idx)
    {
        int ridx = cell_ids[idx];
        for (int kdx=0; kdx<nvars; ++kdx) {
            int rank_idx = ridx + kdx*nlocal;
            int ele_idx = idx + kdx*neles;
            rank_u[rank_idx] = upts[ele_idx];
            rank_du[rank_idx] = rhs[ele_idx];
        }
    }
}

void Exportranku(FlowSys sys, UCFDReal *u)
{
    int n = sys->nvars * sys->nlocal;
    for (int i=0; i<n; ++i)
        u[i] = sys->u[i];
}

void Exportrankdu(FlowSys sys, UCFDReal *du)
{
    int n = sys->nvars * sys->nlocal;
    for (int i=0; i<n; ++i)
        du[i] = sys->du[i];
}

void Exportrankdiag(FlowSys sys, UCFDReal *diag)
{
    int n = sys->nvars * sys->nlocal;
    for (int i=0; i<n; ++i)
        diag[i] = sys->diag[i];
}


/**
 * Pack & update kernels => per-element execution
 */

static ucfd_status_t
UCFDLUSGS_Pack_Impl(const UCFDInt nlocal, const UCFDInt neles, const UCFDInt nvars,
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
    FlowElem *e  = &sys->eles[eidx];

    UCFDCall(UCFDLUSGS_Pack_Impl(
        sys->nlocal, e->neles, sys->nvars, sys->nfvars, turb_factor, a0,
        e->cell_ids, e->uptsb, e->rhs, e->dt, e->dsrc,
        sys->u, sys->du, sys->diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t
UCFDLUSGS_Update_Impl(const UCFDInt nlocal, const UCFDInt neles, const UCFDInt nvars,
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

    UCFDCall(UCFDLUSGS_Update_Impl(
        sys->nlocal, e->neles, sys->nvars,
        e->cell_ids, sys->du, e->uptsb
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}


ucfd_status_t UCFDLUSGS_NSPrepare(FlowSys sys, UCFDReal kappa)
{
    UCFDCall(pre_lusgs(
        sys->nlocal, 0, sys->nfvars, kappa,
        sys->rowptr, sys->slots, sys->face_area, sys->rcp_vol,
        sys->fspr, sys->diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_RANSPrepare(FlowSys sys, UCFDReal kappa)
{
    UCFDCall(pre_lusgs(
        sys->nlocal, sys->nfvars, sys->nvars, kappa,
        sys->rowptr, sys->slots, sys->face_area, sys->rcp_vol,
        sys->tfspr, sys->diag
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_NSLowerSweep(FlowSys sys, UCFDReal kappa)
{
    fluxfunc f = ns_flux_container;
    
    UCFDCall(lower_sweep(
        0, sys->nfvars, kappa, f, sys->fspr,
        sys->nlocal, sys->nvars, sys->nfvars, sys->ndims, sys->nfaces,
        sys->rowptr, sys->colidx, sys->sides, sys->slots, sys->face_area,
        sys->face_normal, sys->rcp_vol, sys->diag,
        sys->u, sys->du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_RANSLowerSweep(FlowSys sys, UCFDReal kappa)
{
#if defined(DEBUG)
    UCFDCheckNull(sys->tfspr, "Turbulent spectral radius is not set\n");
#endif
    fluxfunc f = rans_flux_container;
    
    UCFDCall(lower_sweep(
        sys->nfvars, sys->nvars, kappa, f, sys->tfspr,
        sys->nlocal, sys->nvars, sys->nfvars, sys->ndims, sys->nfaces,
        sys->rowptr, sys->colidx, sys->sides, sys->slots, sys->face_area,
        sys->face_normal, sys->rcp_vol, sys->diag,
        sys->u, sys->du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_NSUpperSweep(FlowSys sys, UCFDReal kappa)
{
    fluxfunc f = ns_flux_container;
    
    UCFDCall(upper_sweep(
        0, sys->nfvars, kappa, f, sys->fspr,
        sys->nlocal, sys->nvars, sys->nfvars, sys->ndims, sys->nfaces,
        sys->rowptr, sys->colidx, sys->sides, sys->slots, sys->face_area,
        sys->face_normal, sys->rcp_vol, sys->diag,
        sys->u, sys->du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_RANSUpperSweep(FlowSys sys, UCFDReal kappa)
{
#if defined(DEBUG)
    UCFDCheckNull(sys->tfspr, "Turbulent spectral radius is not set\n");
#endif
    fluxfunc f = rans_flux_container;
    
    UCFDCall(upper_sweep(
        sys->nfvars, sys->nvars, kappa, f, sys->tfspr,
        sys->nlocal, sys->nvars, sys->nfvars, sys->ndims, sys->nfaces,
        sys->rowptr, sys->colidx, sys->sides, sys->slots, sys->face_area,
        sys->face_normal, sys->rcp_vol, sys->diag,
        sys->u, sys->du
    ));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

