#include <stdio.h>
#include <stdlib.h>
#include "flowsys.h"


ucfd_status_t UCFDFlowSysDestroy(FlowSys *sys)
{
    if (!sys || !*sys) UCFDFunctionReturn(UCFD_SUCCESS);
    UCFDCall((*sys)->destroy(*sys));
    free((*sys)->u);
    free((*sys)->du);
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
                                UCFDReal *rcp_vol)
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
    s->data             = NULL;
    s->destroy          = NULL;

    /* Allocate interior arrays */
    s->u = malloc(nvars*nlocal*sizeof(UCFDReal));
    s->du = malloc(nvars*nlocal*sizeof(UCFDReal));

    /* Elements allocation */
    s->eles = calloc(nelem, sizeof(FlowElem));

    /* Return the newly created system to the caller. */
    *sys = s;

    UCFDFunctionReturn(UCFD_SUCCESS);
}


static ucfd_status_t LUSGSDestroy(FlowSys sys)
{
    LUSGSSys *lusgs = (LUSGSSys *)sys->data;
    free(lusgs->diag);
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDFlowSysSetLUSGS(FlowSys *sys, UCFDReal *fspr, UCFDReal *tfspr)
{
    FlowSys s       = *sys;
    LUSGSSys *lusgs = calloc(1, sizeof(*lusgs));
    UCFDInt nlocal  = s->nlocal;
    UCFDInt nvars   = s->nvars;

    lusgs->diag     = malloc(nvars*nlocal*sizeof(UCFDReal));
    lusgs->fspr     = fspr;
    lusgs->tfspr    = tfspr;

    s->data         = lusgs;
    s->destroy      = LUSGSDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t BLUSGSDestroy(FlowSys sys)
{
    BLUSGSSys *blusgs = (BLUSGSSys *)sys->data;
    free(blusgs->diag);
    free(blusgs->tdiag);
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDFlowSysSetBLUSGS(FlowSys *sys, UCFDReal *jmat, UCFDReal *tjmat)
{
    FlowSys s       = *sys;
    BLUSGSSys *blu  = calloc(1, sizeof(*blu));
    UCFDInt nlocal  = s->nlocal;
    UCFDInt nfvars  = s->nfvars;
    UCFDInt ntvars  = s->nturbvars;

    /* diag : [nlocal, nfvars, nfvars] */
    blu->diag       = malloc(nlocal*nfvars*nfvars*sizeof(UCFDReal));
    if (ntvars != 0)
        blu->tdiag  = malloc(nlocal*ntvars*ntvars*sizeof(UCFDReal));
    else blu->tdiag = NULL;

    blu->jmat       = jmat;
    blu->tjmat      = tjmat;

    s->data         = blu;
    s->destroy      = BLUSGSDestroy;

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
