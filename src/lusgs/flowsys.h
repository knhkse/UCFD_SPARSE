#pragma once

#include "ucfdlusgs.h"
#include "flux.h"


/* Element-level allocation */
typedef struct {
    /* Input from main program */
    UCFDInt neles, nface;

    UCFDReal    *uptsb, *rhs;
    UCFDReal    *dt, *dsrc;
    UCFDInt     *cell_ids;
} FlowElem;


/* Rank-level allocation */
struct _FlowSys {
    /* Basic properties */
    UCFDInt nlocal;
    UCFDInt nvars, nfvars, nturbvars;
    UCFDInt ndims, nfaces;

    /* Sparse matrix format */
    UCFDInt *rowptr;            /* rank_face_indptr */
    UCFDInt *colidx;            /* rank_face_neighbors */
    UCFDInt *slots;             /* rank_face_slots */
    UCFDInt8 *sides;            /* rank_face_sides */
    
    /* Flow attributes */
    UCFDReal *face_area;
    UCFDReal *face_normal;
    UCFDReal *rcp_vol;
    UCFDReal *fspr, *tfspr;

    /* Interior flow attributes */
    UCFDReal *u, *du, *diag;

    /* Elements in own rank */
    FlowElem *eles;
};

UCFD_INTERN ucfd_status_t pre_lusgs(FlowSys, UCFDInt, UCFDInt, UCFDReal, UCFDReal*);
UCFD_INTERN ucfd_status_t lower_sweep(FlowSys, UCFDInt, UCFDInt, UCFDReal, fluxfunc, UCFDReal*);
UCFD_INTERN ucfd_status_t upper_sweep(FlowSys, UCFDInt, UCFDInt, UCFDReal, fluxfunc, UCFDReal*);
