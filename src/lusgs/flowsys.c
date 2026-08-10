#include <stdio.h>
#include <stdlib.h>

#include "flowsys.h"


ucfd_status_t pre_lusgs(FlowSys sys, UCFDInt nv0, UCFDInt nv1, UCFDReal kappa, UCFDReal *lambdaf)
{
    const UCFDInt nlocal = sys->nlocal;
    const UCFDInt *face_indptr = sys->rowptr, *face_slots = sys->slots;
    const UCFDReal *face_area = sys->face_area, *rcp_vol = sys->rcp_vol;
    UCFDReal *diag = sys->diag;

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


static inline void diff_flux(UCFDInt nvars, UCFDInt nfvars, UCFDInt nturbvars, UCFDInt ndims, UCFDInt dnv,
                             fluxfunc fluxf,
                             UCFDReal *u, UCFDReal *du, UCFDReal *df, UCFDReal *nf)
{
    UCFDInt i;
    UCFDReal f[dnv];
    for (i=0; i<nvars; ++i) du[i] += u[i];
    fluxf(nfvars, nturbvars, ndims, u, nf, f);
    fluxf(nfvars, nturbvars, ndims, du, nf, df);
    for (i=0; i<dnv; ++i) df[i] -= f[i];
}


ucfd_status_t lower_sweep(FlowSys sys, UCFDInt nv0, UCFDInt nv1, UCFDReal kappa,
                          fluxfunc fluxf, UCFDReal *lambdaf)
{
    const UCFDInt nvars = sys->nvars, nfvars = sys->nfvars, ntvars = sys->nturbvars;
    const UCFDInt ndims = sys->ndims, nfaces = sys->nfaces, nlocal = sys->nlocal;
    const UCFDInt *face_indptr = sys->rowptr, *face_slots = sys->slots;
    const UCFDInt8 *face_sides = sys->sides;
    const UCFDInt *face_neighbors = sys->colidx;
    const UCFDReal *face_area = sys->face_area, *face_normal = sys->face_normal, *rcp_vol = sys->rcp_vol;
    UCFDReal *rank_u = sys->u, *rank_du = sys->du, *diag = sys->diag;

    UCFDInt ridx, kdx, pos, neib, slot;
    UCFDInt8 side;
    const UCFDInt dnv = nv1 - nv0;
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


ucfd_status_t upper_sweep(FlowSys sys, UCFDInt nv0, UCFDInt nv1, UCFDReal kappa,
                          fluxfunc fluxf, UCFDReal *lambdaf)
{
    const UCFDInt nvars = sys->nvars, nfvars = sys->nfvars, ntvars = sys->nturbvars;
    const UCFDInt ndims = sys->ndims, nfaces = sys->nfaces, nlocal = sys->nlocal;
    const UCFDInt *face_indptr = sys->rowptr, *face_slots = sys->slots;
    const UCFDInt8 *face_sides = sys->sides;
    const UCFDInt *face_neighbors = sys->colidx;
    const UCFDReal *face_area = sys->face_area, *face_normal = sys->face_normal, *rcp_vol = sys->rcp_vol;
    UCFDReal *rank_u = sys->u, *rank_du = sys->du, *diag = sys->diag;

    UCFDInt ridx, kdx, pos, neib, slot;
    UCFDInt8 side;
    const UCFDInt dnv = nv1 - nv0;
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


ucfd_status_t UCFDFlowSysDestroy(FlowSys *sys)
{
    if (!sys || !*sys) UCFDFunctionReturn(UCFD_SUCCESS);
    free((*sys)->u);
    free((*sys)->du);
    free((*sys)->diag);
    free((*sys)->eles);
    free(*sys);
    *sys = NULL;
    
    UCFDFunctionReturn(UCFD_SUCCESS);
}


ucfd_status_t UCFDFlowSysCreate(FlowSys *sys,
                                UCFDInt nlocal, UCFDInt nvars, UCFDInt nfvars,
                                UCFDInt ndims, UCFDInt nfaces, UCFDInt nelem,
                                UCFDInt *rowptr, UCFDInt *colidx,
                                UCFDInt *slots, UCFDInt8 *sides,
                                UCFDReal *face_area, UCFDReal *face_normal,
                                UCFDReal *rcp_vol, UCFDReal *fspr)
{
    FlowSys s = (FlowSys)calloc(1, sizeof(*s));
    UCFDCheckNull(s, "Flow system allocation failed\n");

    /* Get from input arguments */
    s->nlocal           = nlocal;
    s->nvars            = nvars;
    s->nfvars           = nfvars;
    s->nturbvars        = nvars - nfvars;
    s->ndims            = ndims;
    s->nfaces           = nfaces;

    s->rowptr           = rowptr;
    s->colidx           = colidx;
    s->slots            = slots;
    s->sides            = sides;
    s->face_area        = face_area;
    s->face_normal      = face_normal;
    s->rcp_vol          = rcp_vol;
    s->fspr             = fspr;
    s->tfspr            = NULL;

    /* Allocate interior arrays */
    s->u = malloc(nvars*nlocal*sizeof(UCFDReal));
    s->du = malloc(nvars*nlocal*sizeof(UCFDReal));
    s->diag = malloc(nvars*nlocal*sizeof(UCFDReal));

    /* Elements allocation */
    s->eles = calloc(nelem, sizeof(FlowElem));

    /* Return the newly created system to the caller. */
    *sys = s;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDFlowSysSetRANS(FlowSys *sys, UCFDReal *tfspr)
{
    FlowSys s   = *sys;
    s->tfspr    = tfspr;
    UCFDFunctionReturn(UCFD_SUCCESS);
}


ucfd_status_t UCFDFlowSysSetElement(FlowSys *sys, UCFDInt eidx,
                                    UCFDInt neles, UCFDInt nface, UCFDInt *cell_ids,
                                    UCFDReal *uptsb, UCFDReal *rhs,
                                    UCFDReal *dt, UCFDReal *dsrc)
{
    FlowSys s = *sys;
    FlowElem *e = &s->eles[eidx];

    e->neles    = neles;
    e->nface    = nface;
    e->cell_ids = cell_ids;
    e->uptsb    = uptsb;
    e->rhs      = rhs;
    e->dt       = dt;
    e->dsrc     = dsrc;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
