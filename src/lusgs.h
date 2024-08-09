#ifndef LUSGS_H
#define LUSGS_H

/**
 * @file        lusgs.h
 * @brief       Header file for serial LU-SGS method
 * @details     Declaration of each function used in LU-SGS method.
 */

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
void serial_pre_lusgs(int neles, int nface, double factor, \
                      double *fnorm_vol, double *dt, double *diag, double *fspr);


/**
 * @brief       Lower sweep of LU-SGS method for Navier-Stokes equations.
 * @param       neles       Number of element cells
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       ndims       Dimensions
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       mapping     Reordered index by Reverse Cuthill-McKee algorithm [neles]
 * @param       unmapping   Original index before Reverse Cuthill-McKee algorithm [neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       vec_fnorm   Surface vector [nface, ndims, neles]
 * @param       uptsb       Solution array [nfvars, neles]
 * @param       rhsb        Residual (RHS) array [nfvars, neles]
 * @param       dub         Difference array for update [nfvars, neles] (output)
 * @param       diag        Diagonal matrix array [neles]
 * @param       fspr        Wave speed for each cell face [nface, neles]
 */
void ns_serial_lower_sweep(int neles, int nfvars, int nface, int ndims, \
                           int *nei_ele, int *mapping, int *unmapping, double *fnorm_vol, double *vec_fnorm, \
                           double *uptsb, double *rhsb, double *dub, double *diag, double *fspr);


/**
 * @brief       Lower sweep of LU-SGS method for RANS equations.
 * @param       neles       Number of element cells
 * @param       nvars       Number of conservative variables
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       ndims       Dimensions
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       mapping     Reordered index by Reverse Cuthill-McKee algorithm [neles]
 * @param       unmapping   Original index before Reverse Cuthill-McKee algorithm [neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       vec_fnorm   Surface vector [nface, ndims, neles]
 * @param       uptsb       Solution array [nvars, neles]
 * @param       rhsb        Residual (RHS) array [nvars, neles]
 * @param       dub         Difference array for update [nvars, neles] (output)
 * @param       diag        Diagonal matrix array [neles]
 * @param       fspr        Wave speed for each cell face [nface, neles]
 * @param       dsrc        Source term derivatives [nvars, neles]
 */
void rans_serial_lower_sweep(int neles, int nvars, int nfvars, int nface, int ndims, \
                             int *nei_ele, int *mapping, int *unmapping, double *fnorm_vol, double *vec_fnorm, \
                             double *uptsb, double *rhsb, double *dub, double *diag, double *fspr, double *dsrc);


/**
 * @brief       Upper sweep of LU-SGS method for Navier-Stokes equations.
 * @param       neles       Number of element cells
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       ndims       Dimensions
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       mapping     Reordered index by Reverse Cuthill-McKee algorithm [neles]
 * @param       unmapping   Original index before Reverse Cuthill-McKee algorithm [neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       vec_fnorm   Surface vector [nface, ndims, neles]
 * @param       uptsb       Solution array [nfvars, neles]
 * @param       rhsb        Residual (RHS) array [nfvars, neles]
 * @param       dub         Difference array for update [nfvars, neles] (output)
 * @param       diag        Diagonal matrix array [neles]
 * @param       fspr        Wave speed for each cell face [nface, neles]
 */
void ns_serial_upper_sweep(int neles, int nfvars, int nface, int ndims, \
                           int *nei_ele, int *mapping, int *unmapping, double *fnorm_vol, double *vec_fnorm, \
                           double *uptsb, double *rhsb, double *dub, double *diag, double *fspr);


/**
 * @brief       Upper sweep of LU-SGS method for RANS equations.
 * @param       neles       Number of element cells
 * @param       nvars       Number of conservative variables
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       ndims       Dimensions
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       mapping     Reordered index by Reverse Cuthill-McKee algorithm [neles]
 * @param       unmapping   Original index before Reverse Cuthill-McKee algorithm [neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       vec_fnorm   Surface vector [nface, ndims, neles]
 * @param       uptsb       Solution array [nvars, neles]
 * @param       rhsb        Residual (RHS) array [nvars, neles]
 * @param       dub         Difference array for update [nvars, neles] (output)
 * @param       diag        Diagonal matrix array [neles]
 * @param       fspr        Wave speed for each cell face [nface, neles]
 * @param       dsrc        Source term derivatives [nvars, neles]
 */
void rans_serial_upper_sweep(int neles, int nvars, int nfvars, int nface, int ndims, \
                             int *nei_ele, int *mapping, int *unmapping, double *fnorm_vol, double *vec_fnorm, \
                             double *uptsb, double *rhsb, double *dub, double *diag, double *fspr, double *dsrc);


/**
 * @brief       Updates solution array.
 * @param       neles       Number of element cells
 * @param       nvars       Number of conservative variables
 * @param       uptsb       Solution array
 * @param       rhsb        Result of LU-SGS sweeps
 */
void serial_update(int neles, int nvars, double *uptsb, double *rhsb);


#endif // LUSGS_H