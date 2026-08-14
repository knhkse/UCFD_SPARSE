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

UCFD_INTERN
ucfd_status_t pre_lusgs(const UCFDInt nlocal, const UCFDInt nv0, const UCFDInt nv1,
                        const UCFDReal kappa,
                        const UCFDInt *restrict face_indptr, const UCFDInt *restrict face_slots,
                        const UCFDReal *restrict face_area, const UCFDReal *restrict rcp_vol,
                        const UCFDReal *restrict lambdaf, UCFDReal *restrict diag);

UCFD_INTERN
ucfd_status_t lower_sweep(const UCFDInt nv0, const UCFDInt nv1, const UCFDReal kappa,
                          fluxfunc fluxf, const UCFDReal *restrict lambdaf,
                          const UCFDInt nlocal, const UCFDInt nvars, const UCFDInt nfvars,
                          const UCFDInt ndims, const UCFDInt nfaces,
                          const UCFDInt *restrict face_indptr, const UCFDInt *restrict face_neighbors,
                          const UCFDInt8 *restrict face_sides, const UCFDInt *restrict face_slots,
                          const UCFDReal *restrict face_area, const UCFDReal *restrict face_normal,
                          const UCFDReal *restrict rcp_vol, const UCFDReal *restrict diag,
                          const UCFDReal *restrict rank_u, UCFDReal *restrict rank_du);

UCFD_INTERN
ucfd_status_t upper_sweep(const UCFDInt nv0, const UCFDInt nv1, const UCFDReal kappa,
                          fluxfunc fluxf, const UCFDReal *restrict lambdaf,
                          const UCFDInt nlocal, const UCFDInt nvars, const UCFDInt nfvars,
                          const UCFDInt ndims, const UCFDInt nfaces,
                          const UCFDInt *restrict face_indptr, const UCFDInt *restrict face_neighbors,
                          const UCFDInt8 *restrict face_sides, const UCFDInt *restrict face_slots,
                          const UCFDReal *restrict face_area, const UCFDReal *restrict face_normal,
                          const UCFDReal *restrict rcp_vol, const UCFDReal *restrict diag,
                          const UCFDReal *restrict rank_u, UCFDReal *restrict rank_du);

