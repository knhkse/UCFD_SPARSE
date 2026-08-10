#include "flowsys.h"
#include "flux.h"


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
ucfd_status_t UCFDLUSGS_Pack(FlowSys sys, UCFDInt eidx,
                             UCFDReal a0, UCFDReal turb_factor)
{
    FlowElem *e  = &sys->eles[eidx];
    const UCFDInt nvars = sys->nvars, nfvars = sys->nfvars;
    const UCFDInt *cell_ids = e->cell_ids;
    UCFDReal *upts = e->uptsb, *rhs = e->rhs, *dt = e->dt, *dsrc = e->dsrc;
    UCFDReal *rank_u = sys->u, *rank_du = sys->du, *rank_diag = sys->diag;

    const UCFDInt neles = e->neles, nlocal = sys->nlocal;
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

ucfd_status_t UCFDLUSGS_Update(FlowSys sys, UCFDInt eidx)
{
    FlowElem e  = sys->eles[eidx];
    const UCFDInt nvars = sys->nvars, nlocal = sys->nlocal, neles = e.neles;
    const UCFDInt *cell_ids = e.cell_ids;
    UCFDReal *upts = e.uptsb, *rank_du = sys->du;

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


ucfd_status_t UCFDLUSGS_NSPrepare(FlowSys sys, UCFDReal kappa)
{
    UCFDCall(pre_lusgs(sys, 0, sys->nfvars, kappa, sys->fspr));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_RANSPrepare(FlowSys sys, UCFDReal kappa)
{
    UCFDCall(pre_lusgs(sys, sys->nfvars, sys->nvars, kappa, sys->tfspr));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_NSLowerSweep(FlowSys sys, UCFDReal kappa)
{
    fluxfunc f = ns_flux_container;
    
    UCFDCall(lower_sweep(sys, 0, sys->nfvars, kappa, f, sys->fspr));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_RANSLowerSweep(FlowSys sys, UCFDReal kappa)
{
#if defined(DEBUG)
    UCFDCheckNull(sys->tfspr, "Turbulent spectral radius is not set\n");
#endif
    fluxfunc f = rans_flux_container;
    
    UCFDCall(lower_sweep(sys, sys->nfvars, sys->nvars, kappa, f, sys->tfspr));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_NSUpperSweep(FlowSys sys, UCFDReal kappa)
{
    fluxfunc f = ns_flux_container;
    
    UCFDCall(upper_sweep(sys, 0, sys->nfvars, kappa, f, sys->fspr));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDLUSGS_RANSUpperSweep(FlowSys sys, UCFDReal kappa)
{
#if defined(DEBUG)
    UCFDCheckNull(sys->tfspr, "Turbulent spectral radius is not set\n");
#endif
    fluxfunc f = rans_flux_container;
    
    UCFDCall(upper_sweep(sys, sys->nfvars, sys->nvars, kappa, f, sys->tfspr));
    UCFDFunctionReturn(UCFD_SUCCESS);
}

