#include <stdio.h>
#include <stdlib.h>
#include "flowsys.h"


/**
 * Common function
 */
ucfd_status_t UCFDSysCreate(UCFDSys *sys,
                            UCFDInt nlocal, UCFDInt nvars, UCFDInt nfvars,
                            UCFDInt ndims, UCFDInt nfaces, UCFDInt nelem,
                            UCFDInt *face_indptr, UCFDInt *face_slots, UCFDInt *face_neighbors,
                            UCFDInt8 *face_sides, UCFDReal *face_area, UCFDReal *rcp_vol)
{
    UCFDSys s = (UCFDSys)calloc(1, sizeof(*s));
    UCFDCheckNull(s, "Flow system allocation failed\n");

    s->nlocal               = nlocal;
    s->nvars                = nvars;
    s->nfvars               = nfvars;
    s->nturbvars            = nvars - nfvars;
    s->ndims                = ndims;
    s->nfaces               = nfaces;
    s->nelem                = nelem;

    /* Static arrays */
    FlowSys *fsys           = (FlowSys *)calloc(1, sizeof(*fsys));
    FlowEle *eles           = (FlowEle *)calloc(nelem, sizeof(FlowEle));
    fsys->face_indptr       = face_indptr;
    fsys->face_slots        = face_slots;
    fsys->face_neighbors    = face_neighbors;
    fsys->face_sides        = face_sides;
    fsys->face_area         = face_area;
    fsys->rcp_vol           = rcp_vol;
    
    s->eles                 = eles;
    s->flowsys              = fsys;
    *sys                    = s;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDColoredSysCreate(UCFDSys *sys,
                                   UCFDInt nlocal, UCFDInt nvars, UCFDInt nfvars,
                                   UCFDInt ndims, UCFDInt nfaces, UCFDInt nelem)
{
    UCFDSys s = (UCFDSys)calloc(1, sizeof(*s));
    UCFDCheckNull(s, "Flow system allocation failed\n");

    /**
     * !Note!
     * In CUDA BLU-SGS, nlocal = roundup(32, nlocal)
     * to match warp size
     */
    s->nlocal               = nlocal;
    s->nvars                = nvars;
    s->nfvars               = nfvars;
    s->nturbvars            = nvars - nfvars;
    s->ndims                = ndims;
    s->nfaces               = nfaces;
    s->nelem                = nelem;

    PFlowSys *fsys          = (PFlowSys *)calloc(1, sizeof(*fsys));
    PFlowEle *eles          = (PFlowEle *)calloc(nelem, sizeof(PFlowEle));
    s->flowsys              = fsys;
    s->eles                 = eles;
    *sys                    = s;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

/* Only one destroy function */
ucfd_status_t UCFDSysDestroy(UCFDSys *sys)
{
    if (!sys || !*sys) UCFDFunctionReturn(UCFD_SUCCESS);

    free((*sys)->flowsys);
    free((*sys)->eles);
    free(*sys);
    *sys = NULL;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
