/**
 * @file        coloredlusgs.h
 * @brief       Header file for Colored LU-SGS method
 * @details     Declaration of each function used in Colored LU-SGS method.
 *              Parameters are explained here.
 */
#ifndef COLOREDLUSGS_H
#define COLOREDLUSGS_H
#include "config.h"

/**
 * @brief       Computes Diagonal matrix for Colored LU-SGS method.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       factor      Multiplication factor for diffusive margin
 * @param       fnorm_vol   Surface magnitude/cell volume
 * @param       dt          Time step
 * @param       diag        Diagonal matrix for LU-SGS method
 * @param       fspr        Wave speed for each cell face
 */
void parallel_pre_lusgs(UCFDInt neles, UCFDInt nface, UCFDReal factor,
                        UCFDReal *fnorm_vol, UCFDReal *dt, UCFDReal *diag, UCFDReal *fspr);


/**
 * @brief       Lower sweep of Colored LU-SGS method for Navier-Stokes equations.
 * @param       n0          First index for coloring set
 * @param       ne          Last index for coloring set
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells
 * @param       icolor      Color Index of each cell
 * @param       lcolor      Color level of each cell
 * @param       fnorm_vol   Surface magnitude/cell volume
 * @param       vec_fnorm   Surface vector
 * @param       uptsb       Solution array
 * @param       dub         Difference array for update (output)
 * @param       diag        Diagonal matrix array
 * @param       fspr        Wave speed for each cell face
 */
void parallel_ns_lower_sweep(UCFDInt n0, UCFDInt ne, UCFDInt neles, UCFDInt nface,
                             UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                             UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr);


/**
 * @brief       Lower sweep of Colored LU-SGS method for RANS equations.
 * @param       n0          First index for coloring set
 * @param       ne          Last index for coloring set
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells
 * @param       icolor      Color Index of each cell
 * @param       lcolor      Color level of each cell
 * @param       fnorm_vol   Surface magnitude/cell volume
 * @param       vec_fnorm   Surface vector
 * @param       uptsb       Solution array
 * @param       dub         Difference array for update (output)
 * @param       diag        Diagonal matrix array
 * @param       fspr        Wave speed for each cell face
 */
void parallel_rans_lower_sweep(UCFDInt n0, UCFDInt ne, UCFDInt neles, UCFDInt nface,
                               UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                               UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr, UCFDReal *dsrc);


/**
 * @brief       Upper sweep of Colored LU-SGS method for Navier-Stokes equations.
 * @param       n0          First index for coloring set
 * @param       ne          Last index for coloring set
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells
 * @param       icolor      Color Index of each cell
 * @param       lcolor      Color level of each cell
 * @param       fnorm_vol   Surface magnitude/cell volume
 * @param       vec_fnorm   Surface vector
 * @param       uptsb       Solution array
 * @param       dub         Difference array for update (output)
 * @param       diag        Diagonal matrix array
 * @param       fspr        Wave speed for each cell face
 */
void parallel_ns_upper_sweep(UCFDInt n0, UCFDInt ne, UCFDInt neles, UCFDInt nface,
                             UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                             UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr);


/**
 * @brief       Upper sweep of Colored LU-SGS method for RANS equations.
 * @param       n0          First index for coloring set
 * @param       ne          Last index for coloring set
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells
 * @param       icolor      Color Index of each cell
 * @param       lcolor      Color level of each cell
 * @param       fnorm_vol   Surface magnitude/cell volume
 * @param       vec_fnorm   Surface vector
 * @param       uptsb       Solution array
 * @param       dub         Difference array for update (output)
 * @param       diag        Diagonal matrix array
 * @param       fspr        Wave speed for each cell face
 */
void parallel_rans_upper_sweep(UCFDInt n0, UCFDInt ne, UCFDInt neles, UCFDInt nface,
                               UCFDInt *nei_ele, UCFDInt *icolor, UCFDInt *lcolor, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                               UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr, UCFDReal *dsrc);


/**
 * @brief       Updates solution array.
 * @param       neles       Number of element cells
 * @param       uptsb       Solution array
 * @param       dub        Result of LU-SGS sweeps
 */
void parallel_ns_lusgs_update(UCFDInt neles, UCFDReal *uptsb, UCFDReal *dub);


/**
 * @brief       Updates solution array.
 * @param       neles       Number of element cells
 * @param       uptsb       Solution array
 * @param       dub        Result of LU-SGS sweeps
 */
void parallel_lusgs_update(UCFDInt neles, UCFDReal *uptsb, UCFDReal *dub);

#endif