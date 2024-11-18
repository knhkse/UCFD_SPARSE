#ifndef BLUSGS_H
#define BLUSGS_H

/**
 * @file        blusgs.h
 * @brief       Header file for serial LU-SGS method
 * @details     Declaration of each function used in LU-SGS method.
 */

/**
 * @brief       Computes Diagonal matrix for LU-SGS method.
 * @param       neles       Number of element cells
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       factor      Multiplication factor for diffusive margin
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       dt          Time step [neles]
 * @param       diag        (Output) Diagonal matrix for LU-SGS method [nfvars, nfvars, neles]
 * @param       fjmat       Inward block operator include flux Jacobian at each cell face [nfvars, nfvars, nface, neles]
 */
void ns_serial_pre_blusgs(int neles, int nfvars, int nface, double factor, \
                      double *fnorm_vol, double *dt, double *diag, double *fjmat);


/**
 * @brief       Computes Diagonal matrix for Block LU-SGS method for RANS equations.
 * @param       neles       Number of element cells
 * @param       nvars       Number of conservative variables
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       factor      Multiplication factor for diffusive margin
 * @param       betast      Constant \f$\beta^*\f$ for kw-SST turbulent model
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       uptsb       Current solution array [nvars, neles]
 * @param       dt          Time step [neles]
 * @param       tdiag       (Output) Turbulent diagonal matrix array [ntvars, ntvars, neles]
 * @param       tjmat       Inward block operator include flux Jacobian at each cell face [ntvars, ntvars, nface, neles]
 * @param       dsrc        Source term derivatives [nvars, neles]
 */
void rans_serial_pre_blusgs(int neles, int nvars, int nfvars, int nface, double factor, double betast, \
                            double *fnorm_vol, double *uptsb, double *dt, double *tdiag, double *tjmat, double *dsrc);

/**
 * @brief       Lower sweep of Block LU-SGS method for Navier-Stokes equations.
 * @param       neles       Number of element cells
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       mapping     Reordered index by Reverse Cuthill-McKee algorithm [neles]
 * @param       unmapping   Original index before Reverse Cuthill-McKee algorithm [neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       rhsb        Residual (RHS) array [nfvars, neles]
 * @param       dub         (Output) Difference array for update [nvars, neles]
 * @param       diag        Diagonal matrix array [nfvars, nfvars, neles]
 * @param       fjmat       Outward block operator include flux Jacobian at each cell face [nfvars, nfvars, nface, neles]
 */
void ns_serial_block_lower_sweep(int neles, int nfvars, int nface, \
                                 int *nei_ele, int *mapping, int *unmapping, double *fnorm_vol, \
                                 double *rhsb, double *dub, double *diag, double *fjmat);


/**
 * @brief       Lower sweep of Block LU-SGS method for RANS equations.
 * @param       neles       Number of element cells
 * @param       nvars       Number of conservative variables
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       mapping     Reordered index by Reverse Cuthill-McKee algorithm [neles]
 * @param       unmapping   Original index before Reverse Cuthill-McKee algorithm [neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       rhsb        Residual (RHS) array [nvars, neles]
 * @param       dub         (Output) Difference array for update [nvars, neles]
 * @param       tdiag       Turbulent diagonal matrix array [ntvars, ntvars, neles]
 * @param       tjmat       Outward block operator include flux Jacobian at each cell face [ntvars, ntvars, nface, neles]
 */
void rans_serial_block_lower_sweep(int neles, int nvars, int nfvars, int nface, \
                                   int *nei_ele, int *mapping, int *unmapping, double *fnorm_vol, \
                                   double *rhsb, double *dub, double *tdiag, double *tjmat);


/**
 * @brief       Upper sweep of Block LU-SGS method for Navier-Stokes equations.
 * @param       neles       Number of element cells
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       mapping     Reordered index by Reverse Cuthill-McKee algorithm [neles]
 * @param       unmapping   Original index before Reverse Cuthill-McKee algorithm [neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       rhsb        Residual (RHS) array [nfvars, neles]
 * @param       dub         (Output) Difference array for update [nvars, neles]
 * @param       diag        Diagonal matrix array [nfvars, nfvars, neles]
 * @param       fjmat       Outward block operator include flux Jacobian at each cell face [nfvars, nfvars, nface, neles]
 */
void ns_serial_block_upper_sweep(int neles, int nfvars, int nface, \
                                 int *nei_ele, int *mapping, int *unmapping, double *fnorm_vol, \
                                 double *rhsb, double *dub, double *diag, double *fjmat);


/**
 * @brief       Upper sweep of Block LU-SGS method for RANS equations.
 * @param       neles       Number of element cells
 * @param       nvars       Number of conservative variables
 * @param       nfvars      Number of flux variables
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       mapping     Reordered index by Reverse Cuthill-McKee algorithm [neles]
 * @param       unmapping   Original index before Reverse Cuthill-McKee algorithm [neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       rhsb        Residual (RHS) array [nvars, neles]
 * @param       dub         (Output) Difference array for update [nvars, neles]
 * @param       tdiag       Turbulent diagonal matrix array [ntvars, ntvars, neles]
 * @param       tjmat       Outward block operator include flux Jacobian at each cell face [ntvars, ntvars, nface, neles]
 */
void rans_serial_block_upper_sweep(int neles, int nvars, int nfvars, int nface, \
                                 int *nei_ele, int *mapping, int *unmapping, double *fnorm_vol, \
                                 double *rhsb, double *dub, double *tdiag, double *tjmat);


/**
 * @brief       Updates solution array.
 * @param       neles       Number of element cells
 * @param       nvars       Number of conservative variables
 * @param       uptsb       Solution array
 * @param       dub         Result of Block LU-SGS sweeps
 * @param       subres      Residual of each sub-iteration
 */
void serial_update(int neles, int nvars, double *uptsb, double *dub, double *subres);


#endif // BLUSGS_H