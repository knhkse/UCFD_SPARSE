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


/* LU-SGS system */
typedef struct {
    UCFDReal *diag;
    UCFDReal *fspr, *tfspr;
} LUSGSSys;

/* Block LU-SGS system */
typedef struct {
    UCFDReal *diag, *tdiag;
    UCFDReal *jmat, *tjmat;
} BLUSGSSys;

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

    /* Interior flow attributes */
    UCFDReal *u, *du;

    /* LU-SGS or BLU-SGS */
    void *data;

    /* Elements in own rank */
    FlowElem *eles;

    /* Destruction function */
    ucfd_status_t (*destroy)(FlowSys);
};
