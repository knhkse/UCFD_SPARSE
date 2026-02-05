/**
 * @file        coloredblusgs.h
 * @brief       Header file for colored Block LU-SGS method
 */
#ifndef COLOREDBLUSGS_H
#define COLOREDBLUSGS_H
#include "config.h"

void ns_parallel_pre_blusgs(UCFD_INT neles, UCFD_INT nface, UCFD_FLOAT factor,
                          UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *dt, UCFD_FLOAT *diag, UCFD_FLOAT *fjmat);


void rans_parallel_pre_blusgs(UCFD_INT neles, UCFD_INT nface, UCFD_FLOAT factor,
                            UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *uptsb, UCFD_FLOAT *dt,
                            UCFD_FLOAT *tdiag, UCFD_FLOAT *tjmat, UCFD_FLOAT *dsrc);


void ns_parallel_block_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                             UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol,
                             UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fjmat);


void rans_parallel_block_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                               UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol,
                               UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *tdiag, UCFD_FLOAT *tjmat);


void blusgs_parallel_update(UCFD_INT neles, UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *subres);


#endif