/**
 * @file        coloredblusgs.h
 * @brief       Header file for colored Block LU-SGS method
 */
#ifndef COLOREDBLUSGS_H
#define COLOREDBLUSGS_H
#include "config.h"

void parallel_ns_pre_blusgs(UCFDInt neles, UCFDInt nface, UCFDReal factor,
                          UCFDReal *fnorm_vol, UCFDReal *dt, UCFDReal *diag, UCFDReal *fjmat);


void parallel_rans_pre_blusgs(UCFDInt neles, UCFDInt nface, UCFDReal factor,
                            UCFDReal *fnorm_vol, UCFDReal *uptsb, UCFDReal *dt,
                            UCFDReal *tdiag, UCFDReal *tjmat, UCFDReal *dsrc);


void parallel_ns_block_sweep(UCFDInt n0, UCFDInt ne, UCFDInt neles, UCFDInt nface,
                             UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol,
                             UCFDReal *rhsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fjmat);


void parallel_rans_block_sweep(UCFDInt n0, UCFDInt ne, UCFDInt neles, UCFDInt nface,
                               UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol,
                               UCFDReal *rhsb, UCFDReal *dub, UCFDReal *tdiag, UCFDReal *tjmat);


void parallel_blusgs_update(UCFDInt neles, UCFDReal *uptsb, UCFDReal *dub);

#endif