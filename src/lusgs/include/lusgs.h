/**
 * @file        lusgs.h
 * @brief       Header file for serial LU-SGS method
 * @details     Declaration of each function used in LU-SGS method.
 */
#ifndef LUSGS_H
#define LUSGS_H
#include "config.h"

/**
 * @brief       Computes Diagonal matrix for LU-SGS method.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       factor      Multiplication factor for diffusive margin
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       dt          Time step [neles]
 * @param       diag        Diagonal matrix for LU-SGS method [neles]
 * @param       fspr        Wave speed for each cell face [nface, neles]
 */
void serial_pre_lusgs(UCFDInt neles, UCFDInt nface, UCFDReal factor,
                      UCFDReal *fnorm_vol, UCFDReal *dt, UCFDReal *diag, UCFDReal *fspr);


/**
 * @brief       Lower sweep of LU-SGS method for Navier-Stokes equations.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       vec_fnorm   Surface vector [nface, ndims, neles]
 * @param       uptsb       Solution array [nfvars, neles]
 * @param       rhsb        Residual (RHS) array [nfvars, neles]
 * @param       dub         Difference array for update [nfvars, neles] (output)
 * @param       diag        Diagonal matrix array [neles]
 * @param       fspr        Wave speed for each cell face [nface, neles]
 */
void serial_ns_lower_sweep(UCFDInt neles, UCFDInt nface, UCFDInt *nei_ele, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                           UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr);


/**
 * @brief       Lower sweep of LU-SGS method for RANS equations.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       vec_fnorm   Surface vector [nface, ndims, neles]
 * @param       uptsb       Solution array [nvars, neles]
 * @param       dub         Difference array for update [nvars, neles] (output)
 * @param       diag        Diagonal matrix array [neles]
 * @param       fspr        Wave speed for each cell face [nface, neles]
 * @param       dsrc        Source term derivatives [nvars, neles]
 */
void serial_rans_lower_sweep(UCFDInt neles, UCFDInt nface, UCFDInt *nei_ele, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm, 
                             UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr, UCFDReal *dsrc);


/**
 * @brief       Upper sweep of LU-SGS method for Navier-Stokes equations.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       vec_fnorm   Surface vector [nface, ndims, neles]
 * @param       uptsb       Solution array [nfvars, neles]
 * @param       rhsb        Residual (RHS) array [nfvars, neles]
 * @param       dub         Difference array for update [nfvars, neles] (output)
 * @param       diag        Diagonal matrix array [neles]
 * @param       fspr        Wave speed for each cell face [nface, neles]
 */
void serial_ns_upper_sweep(UCFDInt neles, UCFDInt nface, UCFDInt *nei_ele, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                           UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr);


/**
 * @brief       Upper sweep of LU-SGS method for RANS equations.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       vec_fnorm   Surface vector [nface, ndims, neles]
 * @param       uptsb       Solution array [nvars, neles]
 * @param       dub         Difference array for update [nvars, neles] (output)
 * @param       diag        Diagonal matrix array [neles]
 * @param       fspr        Wave speed for each cell face [nface, neles]
 * @param       dsrc        Source term derivatives [nvars, neles]
 */
void serial_rans_upper_sweep(UCFDInt neles, UCFDInt nface, UCFDInt *nei_ele, UCFDReal *fnorm_vol, UCFDReal *vec_fnorm,
                             UCFDReal *uptsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fspr, UCFDReal *dsrc);


/**
 * @brief       Updates solution array.
 * @param       neles       Number of element cells
 * @param       uptsb       Solution array
 * @param       rhsb        Result of LU-SGS sweeps
 */
void serial_lusgs_update(UCFDInt neles, UCFDReal *uptsb, UCFDReal *rhsb);

#endif // LUSGS_H