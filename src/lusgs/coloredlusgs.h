#ifndef COLOREDLUSGS_H
#define COLOREDLUSGS_H

/**
 * @file        coloredlusgs.h
 * @brief       Header file for Colored LU-SGS method
 * @details     Declaration of each function used in Colored LU-SGS method.
 *              Parameters are explained here.
 */

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
void parallel_pre_lusgs(int neles, int nface, double factor, \
                        double *fnorm_vol, double *dt, double *diag, double *fspr);


/**
 * @brief       Lower sweep of Colored LU-SGS method for Navier-Stokes equations.
 * @param       n0          First index for coloring set
 * @param       ne          Last index for coloring set
 * @param       neles       Number of element cells
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       ndims       Dimensions
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
void ns_parallel_lower_sweep(int n0, int ne, int neles, int nfvars, int nface, int ndims, \
                             int *nei_ele, int *icolor, int *lcolor, double *fnorm_vol, double *vec_fnorm, \
                             double *uptsb, double *rhsb, double *dub, double *diag, double *fspr);


/**
 * @brief       Lower sweep of Colored LU-SGS method for RANS equations.
 * @param       n0          First index for coloring set
 * @param       ne          Last index for coloring set
 * @param       neles       Number of element cells
 * @param       nvars       Number of conservative variables
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       ndims       Dimensions
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
void rans_parallel_lower_sweep(int n0, int ne, int neles, int nvars, int nfvars, int nface, int ndims, \
                               int *nei_ele, int *icolor, int *lcolor, double *fnorm_vol, double *vec_fnorm, \
                               double *uptsb, double *rhsb, double *dub, double *diag, double *fspr, double *dsrc);


/**
 * @brief       Upper sweep of Colored LU-SGS method for Navier-Stokes equations.
 * @param       n0          First index for coloring set
 * @param       ne          Last index for coloring set
 * @param       neles       Number of element cells
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       ndims       Dimensions
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
void ns_parallel_upper_sweep(int n0, int ne, int neles, int nfvars, int nface, int ndims, \
                             int *nei_ele, int *icolor, int *lcolor, double *fnorm_vol, double *vec_fnorm, \
                             double *uptsb, double *rhsb, double *dub, double *diag, double *fspr);


/**
 * @brief       Upper sweep of Colored LU-SGS method for RANS equations.
 * @param       n0          First index for coloring set
 * @param       ne          Last index for coloring set
 * @param       neles       Number of element cells
 * @param       nvars       Number of conservative variables
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       ndims       Dimensions
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
void rans_parallel_upper_sweep(int n0, int ne, int neles, int nvars, int nfvars, int nface, int ndims, \
                               int *nei_ele, int *icolor, int *lcolor, double *fnorm_vol, double *vec_fnorm, \
                               double *uptsb, double *rhsb, double *dub, double *diag, double *fspr, double *dsrc);


/**
 * @brief       Updates solution array.
 * @param       neles       Number of element cells
 * @param       nvars       Number of conservative variables
 * @param       uptsb       Solution array
 * @param       rhsb        Result of LU-SGS sweeps
 */
void parallel_update(int neles, int nvars, double *uptsb, double *rhsb);


#endif