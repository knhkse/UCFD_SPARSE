#pragma once

#include "ucfdlusgs.h"
#include "flux.h"


/* LU-SGS system */
typedef struct {
    UCFDReal *u, *du;
    UCFDReal *diag;
    UCFDReal *fspr, *tfspr;
} LUSGSSys;

/* Block LU-SGS system */
typedef struct {
    UCFDReal *rhs, *du, *dup;
    UCFDReal *diag, *tdiag;
    UCFDReal *jmat, *tjmat;
} BLUSGSSys;

/* CUDA BLU-SGS system */
typedef struct {
    UCFDReal *du, *dup;
    UCFDReal *diag, *tdiag;
    UCFDReal *jmat, *tjmat;
} CUBLUSGSSys;


/**
 * Serial execution system
 */
typedef struct {
    /* Input from main program */
    UCFDInt     neles, nface;

    UCFDInt     *cell_ids;      /* rank_cell_ids */
    UCFDReal    *uptsb, *rhs;
    UCFDReal    *dt, *dsrc;
    UCFDReal    *vol, *resid_out;
} FlowElem;

struct _FlowSys {
    /* Basic properties */
    UCFDInt     nlocal;
    UCFDInt     nvars, nfvars, nturbvars;
    UCFDInt     ndims, nfaces;

    /* Sparse matrix format */
    UCFDInt     *rowptr;            /* rank_face_indptr */
    UCFDInt     *colidx;            /* rank_face_neighbors */
    UCFDInt     *slots;             /* rank_face_slots */
    UCFDInt8    *sides;            /* rank_face_sides */
    
    /* Flow attributes */
    UCFDReal    *face_area;
    UCFDReal    *face_normal;
    UCFDReal    *rcp_vol;

    /* Specific relaxation */
    void        *data;

    /* Elements in own rank */
    FlowElem    *eles;

    /* Destruction function */
    ucfd_status_t (*destroy)(FlowSys);
};


/**
 * Parallel execution system
 * use rank-coloring grid
 */
typedef struct {
    FlowElem    base;           /* rank_ele_cell_ids */

    /* Coloring attributes */
    UCFDInt     ncolors, *nblocks;
    UCFDInt     *face_refs;
    UCFDReal    *face_factors;
    UCFDInt     *color_offsets, *color_order;

    /* Colored LU-SGS */
    UCFDInt     *lower_neighbors, *upper_neighbors;
    UCFDReal    *face_normal;
    UCFDReal    *wave_factors;      /* Internal allocation */

    /* Colored BLU-SGS */
    UCFDInt     *neighbors;
    /* offdiag : [neles, nface, nfvars, nfvars] */
    /* turb_offdiag : [neles, nface, nturbvars, nturbvars] */
    UCFDReal    *flow_offdiag, *turb_offdiag;   /* Internal allocation */
    
} PFlowElem;

struct _PFlowSys {
    /* Basic properties */
    UCFDInt     nlocal, nelem;
    UCFDInt     nvars, nfvars, nturbvars;
    UCFDInt     ndims, nfaces;

    /* LU-SGS or BLU-SGS */
    void        *data;

    /* Elements in own rank */
    PFlowElem   *eles;

    /* Destruction function */
    ucfd_status_t (*destroy)(PFlowSys);
};
