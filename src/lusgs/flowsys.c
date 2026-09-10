#include <stdio.h>
#include <stdlib.h>
#include "flowsys.h"


/**
 * rank-order mesh
 */
ucfd_status_t UCFDPFlowSysDestroy(PFlowSys *sys)
{
    if (!sys || !*sys) UCFDFunctionReturn(UCFD_SUCCESS);
    UCFDCall((*sys)->destroy(*sys));

    free((*sys)->data);
    free((*sys)->eles);
    free(*sys);
    *sys = NULL;
    
    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDFlowSysDestroy(FlowSys *sys)
{
    if (!sys || !*sys) UCFDFunctionReturn(UCFD_SUCCESS);
    UCFDCall((*sys)->destroy(*sys));

    free((*sys)->data);
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

    /* Elements allocation */
    s->eles = (FlowElem *)calloc(nelem, sizeof(FlowElem));

    /* Return the newly created system to the caller. */
    *sys = s;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDFlowSysSetElement(FlowSys sys, UCFDInt eidx,
                                    UCFDInt neles, UCFDInt nface, UCFDInt *cell_ids,
                                    UCFDReal *uptsb, UCFDReal *rhs,
                                    UCFDReal *dt, UCFDReal *dsrc,
                                    UCFDReal *vol, UCFDReal *resid_out)
{
    FlowElem *eles = sys->eles;
    FlowElem *e = &eles[eidx];

    e->neles        = neles;
    e->nface        = nface;
    e->cell_ids     = cell_ids;
    e->uptsb        = uptsb;
    e->rhs          = rhs;
    e->dt           = dt;
    e->dsrc         = dsrc;
    e->vol          = vol;
    e->resid_out    = resid_out;

    UCFDFunctionReturn(UCFD_SUCCESS);
}


/**
 * rank-coloring mesh
 */
ucfd_status_t UCFDPFlowSysCreate(PFlowSys *sys,
                                 UCFDInt nlocal, UCFDInt nvars, UCFDInt nfvars,
                                 UCFDInt ndims, UCFDInt nfaces, UCFDInt nelem)
{
    PFlowSys s = (PFlowSys)calloc(1, sizeof(*s));
    UCFDCheckNull(s, "Parallel Flow system allocation failed\n");

    /* Get from input arguments */
    s->nlocal           = nlocal;
    s->nelem            = nelem;
    s->nvars            = nvars;
    s->nfvars           = nfvars;
    s->nturbvars        = nvars - nfvars;
    s->ndims            = ndims;
    s->nfaces           = nfaces;
    s->data             = NULL;
    s->destroy          = NULL;

    /* Element allocation */
    s->eles = calloc(nelem, sizeof(PFlowElem));

    /* Return the newly created system to the caller. */
    *sys = s;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
