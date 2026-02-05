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
void parallel_pre_lusgs(UCFD_INT neles, UCFD_INT nface, UCFD_FLOAT factor,
                        UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *dt, UCFD_FLOAT *diag, UCFD_FLOAT *fspr);


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
 * @param       rhsb        Residual (RHS) array
 * @param       dub         Difference array for update (output)
 * @param       diag        Diagonal matrix array
 * @param       fspr        Wave speed for each cell face
 */
void ns_parallel_lower_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                             UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                             UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr);


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
 * @param       rhsb        Residual (RHS) array
 * @param       dub         Difference array for update (output)
 * @param       diag        Diagonal matrix array
 * @param       fspr        Wave speed for each cell face
 */
void rans_parallel_lower_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                               UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                               UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr, UCFD_FLOAT *dsrc);


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
 * @param       rhsb        Residual (RHS) array
 * @param       dub         Difference array for update (output)
 * @param       diag        Diagonal matrix array
 * @param       fspr        Wave speed for each cell face
 */
void ns_parallel_upper_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                             UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                             UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr);


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
 * @param       rhsb        Residual (RHS) array
 * @param       dub         Difference array for update (output)
 * @param       diag        Diagonal matrix array
 * @param       fspr        Wave speed for each cell face
 */
void rans_parallel_upper_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                               UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                               UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr, UCFD_FLOAT *dsrc);


/**
 * @brief       Updates solution array.
 * @param       neles       Number of element cells
 * @param       uptsb       Solution array
 * @param       rhsb        Result of LU-SGS sweeps
 */
void clusgs_parallel_update(UCFD_INT neles, UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb);

#endif