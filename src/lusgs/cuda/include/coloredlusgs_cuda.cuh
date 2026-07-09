#ifndef COLOREDLUSGS_CUDA_H
#define COLOREDLUSGS_CUDA_H
#include "config.h"


__global__ void
pre_lusgs_kernel(UCFDInt neles, UCFDInt nface, UCFDReal factor,
                 UCFDReal *fnorm_vol, UCFDReal *dt, UCFDReal *diag, UCFDReal *fspr);


__global__ void
ns_lower_sweep_kernel(UCFDInt interval, UCFDInt n0, UCFDInt neles, UCFDInt nface,
                      UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                      UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr);


__global__ void
rans_lower_sweep_kernel(UCFDInt interval, UCFDInt n0, UCFDInt neles, UCFDInt nface,
                        UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                        UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr, UCFDReal *dsrc);


__global__ void
ns_upper_sweep_kernel(UCFDInt interval, UCFDInt n0, UCFDInt neles, UCFDInt nface,
                      UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                      UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr);


__global__ void
rans_upper_sweep_kernel(UCFDInt interval, UCFDInt n0, UCFDInt neles, UCFDInt nface,
                        UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                        UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr, UCFDReal *dsrc);


__global__ void
lusgs_update_kernel(UCFDInt neles, UCFDReal *uptsb, UCFDReal *dub);


extern "C" void
cuda_pre_lusgs(UCFDInt neles, UCFDInt nface, UCFDReal factor,
               UCFDReal *fnorm_vol, UCFDReal *dt, UCFDReal *diag, UCFDReal *fspr);

extern "C" void
cuda_ns_lower_sweep(UCFDInt n0, UCFDInt ne, UCFDInt neles, UCFDInt nface,
                    UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                    UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr);

extern "C" void
cuda_rans_lower_sweep(UCFDInt n0, UCFDInt ne, UCFDInt neles, UCFDInt nface,
                      UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                      UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr, UCFDReal *dsrc);

extern "C" void
cuda_ns_upper_sweep(UCFDInt n0, UCFDInt ne, UCFDInt neles, UCFDInt nface,
                    UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                    UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr);

extern "C" void
cuda_rans_upper_sweep(UCFDInt n0, UCFDInt ne, UCFDInt neles, UCFDInt nface,
                      UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                      UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr, UCFDReal *dsrc);

extern "C" void
cuda_lusgs_update(UCFDInt neles, UCFDReal *uptsb, UCFDReal *dub);

#endif