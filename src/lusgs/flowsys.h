#pragma once

#include "ucfdlusgs.h"
#include "flux.h"

/**
 * Entire system
 */
struct _UCFDSys {
    /* Basic scalar properties */
    UCFDInt     nlocal, nelem;
    UCFDInt     nvars, nfvars, nturbvars;
    UCFDInt     ndims, nfaces;

    /* Specific system structure
        FlowSys or PFlowSys
    */
    void        *flowsys;
    void        *eles;
};

/**
 * Serial / parallel system
 */
#define BASESYS                                 \
    /**                                         \
     * Basic system arrays                      \
     * LU-SGS : rank_u, rank_du, rank_diag      \
     * BLU-SGS : rank_du, rank_dup, rank_rhs    \
     */                                         \
    UCFDReal    *rank_u, *rank_du, *rank_dup;   \
    UCFDReal    *rank_diag, *rank_rhs;          \
\
    /* LU-SGS */                                \
    UCFDReal    *fspr, *tfspr;                  \
\
    /* BLU-SGS */                               \
    UCFDReal    *flow_diag, *turb_diag;         \
    UCFDReal    *jmat, *tjmat;

typedef struct {
    BASESYS
    UCFDInt     *face_indptr, *face_slots, *face_neighbors;
    UCFDInt8    *face_sides;
    UCFDReal    *face_area, *face_normal, *rcp_vol;
} FlowSys;

typedef struct {
    BASESYS
} PFlowSys;




/**
 * Note :
 * Colored element has loop-contained arrays
 */
#define BASEELEMENT                             \
    UCFDInt     neles, nface;                   \
    UCFDInt     *cell_ids;                      \
    UCFDReal    *uptsb, *rhs, *resid_out;       \
    UCFDReal    *dt, *dsrc, *vol;

typedef struct {
    BASEELEMENT
} FlowEle;

typedef struct {
    BASEELEMENT
    UCFDInt     *face_refs, *neighbors;
    UCFDReal    *face_factors;
    UCFDReal    *flow_offdiag, *turb_offdiag;

    /* Coloring properties */
    UCFDInt     cell_offset;
    UCFDInt     ncolors, *nblocks;
    UCFDInt     *color_offset;
} PFlowEle;


/* Internal functions */
UCFD_INTERN ucfd_status_t _rank_blusgs_update(const UCFDInt, const UCFDInt, const UCFDInt,
                                              const UCFDInt*, const UCFDReal*, UCFDReal*);
UCFD_INTERN ucfd_status_t _rank_tblusgs_pack(const UCFDInt, const UCFDInt, const UCFDInt,
                                            const UCFDInt, const UCFDReal, const UCFDReal,
                                            srcjacobian, const UCFDInt*, const UCFDReal*,
                                            const UCFDReal*, const UCFDReal*, UCFDReal*);
UCFD_INTERN ucfd_status_t _rank_sub_residual(const UCFDInt, const UCFDInt, const UCFDInt,
                                             const UCFDInt*, const UCFDReal*, const UCFDReal*,
                                             UCFDReal*, UCFDReal*);
