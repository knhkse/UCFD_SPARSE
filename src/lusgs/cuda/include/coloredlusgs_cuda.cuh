#ifndef COLOREDLUSGS_CUDA_H
#define COLOREDLUSGS_CUDA_H
#include "config.h"


__global__ void
pre_lusgs_kernel(UCFD_INT neles, UCFD_INT nface, UCFD_FLOAT factor,
                 UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *dt, UCFD_FLOAT *diag, UCFD_FLOAT *fspr);


__global__ void
ns_lower_sweep_kernel(UCFD_INT interval, UCFD_INT n0, UCFD_INT neles, UCFD_INT nface,
                      UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                      UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr);


__global__ void
rans_lower_sweep_kernel(UCFD_INT interval, UCFD_INT n0, UCFD_INT neles, UCFD_INT nface,
                        UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                        UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr, UCFD_FLOAT *dsrc);


__global__ void
ns_upper_sweep_kernel(UCFD_INT interval, UCFD_INT n0, UCFD_INT neles, UCFD_INT nface,
                      UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                      UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr);


__global__ void
rans_upper_sweep_kernel(UCFD_INT interval, UCFD_INT n0, UCFD_INT neles, UCFD_INT nface,
                        UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                        UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr, UCFD_FLOAT *dsrc);


__global__ void
lusgs_update_kernel(UCFD_INT neles, UCFD_FLOAT *uptsb, UCFD_FLOAT *dub);


extern "C" void
cuda_pre_lusgs(UCFD_INT neles, UCFD_INT nface, UCFD_FLOAT factor,
               UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *dt, UCFD_FLOAT *diag, UCFD_FLOAT *fspr);

extern "C" void
cuda_ns_lower_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                    UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                    UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr);

extern "C" void
cuda_rans_lower_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                      UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                      UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr, UCFD_FLOAT *dsrc);

extern "C" void
cuda_ns_upper_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                    UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                    UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr);

extern "C" void
cuda_rans_upper_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                      UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                      UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr, UCFD_FLOAT *dsrc);

extern "C" void
cuda_lusgs_update(UCFD_INT neles, UCFD_FLOAT *uptsb, UCFD_FLOAT *dub);

#endif