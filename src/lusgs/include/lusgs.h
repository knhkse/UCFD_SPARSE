/**
 * @file        lusgs.h
 * @brief       Header file for serial LU-SGS method
 * @details     Declaration of each function used in LU-SGS method.
 */
#ifndef LUSGS_H
#define LUSGS_H
#include "flux.h"


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
void serial_pre_lusgs(UCFD_INT neles, UCFD_INT nface, UCFD_FLOAT factor,
                      UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *dt, UCFD_FLOAT *diag, UCFD_FLOAT *fspr);


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
void ns_serial_lower_sweep(UCFD_INT neles, UCFD_INT nface, UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                           UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr);


/**
 * @brief       Lower sweep of LU-SGS method for RANS equations.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       vec_fnorm   Surface vector [nface, ndims, neles]
 * @param       uptsb       Solution array [nvars, neles]
 * @param       rhsb        Residual (RHS) array [nvars, neles]
 * @param       dub         Difference array for update [nvars, neles] (output)
 * @param       diag        Diagonal matrix array [neles]
 * @param       fspr        Wave speed for each cell face [nface, neles]
 * @param       dsrc        Source term derivatives [nvars, neles]
 */
void rans_serial_lower_sweep(UCFD_INT neles, UCFD_INT nface, UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm, 
                             UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr, UCFD_FLOAT *dsrc);


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
void ns_serial_upper_sweep(UCFD_INT neles, UCFD_INT nface, UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                           UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr);


/**
 * @brief       Upper sweep of LU-SGS method for RANS equations.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       vec_fnorm   Surface vector [nface, ndims, neles]
 * @param       uptsb       Solution array [nvars, neles]
 * @param       rhsb        Residual (RHS) array [nvars, neles]
 * @param       dub         Difference array for update [nvars, neles] (output)
 * @param       diag        Diagonal matrix array [neles]
 * @param       fspr        Wave speed for each cell face [nface, neles]
 * @param       dsrc        Source term derivatives [nvars, neles]
 */
void rans_serial_upper_sweep(UCFD_INT neles, UCFD_INT nface, UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                             UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr, UCFD_FLOAT *dsrc);


/**
 * @brief       Updates solution array.
 * @param       neles       Number of element cells
 * @param       uptsb       Solution array
 * @param       rhsb        Result of LU-SGS sweeps
 */
void lusgs_serial_update(UCFD_INT neles, UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb);

#endif // LUSGS_H