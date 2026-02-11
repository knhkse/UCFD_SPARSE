/**
 * @file        blusgs.h
 * @brief       Header file for serial LU-SGS method
 * @details     Declaration of each function used in LU-SGS method.
 */
#ifndef BLUSGS_H
#define BLUSGS_H
#include "config.h"

/**
 * @brief       Computes Diagonal matrix for LU-SGS method.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       factor      Multiplication factor for diffusive margin
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       dt          Time step [neles]
 * @param       diag        (Output) Diagonal matrix for LU-SGS method [nfvars, nfvars, neles]
 * @param       fjmat       Inward block operator include flux Jacobian at each cell face [nfvars, nfvars, nface, neles]
 */
void ns_serial_pre_blusgs(UCFD_INT neles, UCFD_INT nface, UCFD_FLOAT factor,
                          UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *dt, UCFD_FLOAT *diag, UCFD_FLOAT *fjmat);


/**
 * @brief       Computes Diagonal matrix for Block LU-SGS method for RANS equations.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       factor      Multiplication factor for diffusive margin
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       uptsb       Current solution array [nvars, neles]
 * @param       dt          Time step [neles]
 * @param       tdiag       (Output) Turbulent diagonal matrix array [ntvars, ntvars, neles]
 * @param       tjmat       Inward block operator include flux Jacobian at each cell face [ntvars, ntvars, nface, neles]
 * @param       dsrc        Source term derivatives [nvars, neles]
 */
void rans_serial_pre_blusgs(UCFD_INT neles, UCFD_INT nface, UCFD_FLOAT factor,
                            UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *uptsb, UCFD_FLOAT *dt,
                            UCFD_FLOAT *tdiag, UCFD_FLOAT *tjmat, UCFD_FLOAT *dsrc);


/**
 * @brief       Lower sweep of Block LU-SGS method for Navier-Stokes equations.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       rhsb        Residual (RHS) array [nfvars, neles]
 * @param       dub         (Output) Difference array for update [nvars, neles]
 * @param       diag        Diagonal matrix array [nfvars, nfvars, neles]
 * @param       fjmat       Outward block operator include flux Jacobian at each cell face [nfvars, nfvars, nface, neles]
 */
void ns_serial_block_lower_sweep(UCFD_INT neles, UCFD_INT nface,
                                 UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol,
                                 UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fjmat);


/**
 * @brief       Lower sweep of Block LU-SGS method for RANS equations.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       rhsb        Residual (RHS) array [nvars, neles]
 * @param       dub         (Output) Difference array for update [nvars, neles]
 * @param       tdiag       Turbulent diagonal matrix array [ntvars, ntvars, neles]
 * @param       tjmat       Outward block operator include flux Jacobian at each cell face [ntvars, ntvars, nface, neles]
 */
void rans_serial_block_lower_sweep(UCFD_INT neles, UCFD_INT nface,
                                   UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol,
                                   UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *tdiag, UCFD_FLOAT *tjmat);


/**
 * @brief       Upper sweep of Block LU-SGS method for Navier-Stokes equations.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       rhsb        Residual (RHS) array [nfvars, neles]
 * @param       dub         (Output) Difference array for update [nvars, neles]
 * @param       diag        Diagonal matrix array [nfvars, nfvars, neles]
 * @param       fjmat       Outward block operator include flux Jacobian at each cell face [nfvars, nfvars, nface, neles]
 */
void ns_serial_block_upper_sweep(UCFD_INT neles, UCFD_INT nface,
                                 UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol,
                                 UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fjmat);


/**
 * @brief       Upper sweep of Block LU-SGS method for RANS equations.
 * @param       neles       Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       fnorm_vol   Surface magnitude/cell volume [nface, neles]
 * @param       rhsb        Residual (RHS) array [nvars, neles]
 * @param       dub         (Output) Difference array for update [nvars, neles]
 * @param       tdiag       Turbulent diagonal matrix array [ntvars, ntvars, neles]
 * @param       tjmat       Outward block operator include flux Jacobian at each cell face [ntvars, ntvars, nface, neles]
 */
void rans_serial_block_upper_sweep(UCFD_INT neles, UCFD_INT nface,
                                   UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol,
                                   UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *tdiag, UCFD_FLOAT *tjmat);


/**
 * @brief       Updates solution array.
 * @param       neles       Number of element cells
 * @param       uptsb       Solution array
 * @param       dub         Result of Block LU-SGS sweeps
 * @param       subres      Residual of each sub-iteration
 */
void blusgs_serial_ns_update(UCFD_INT neles, UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *subres);


/**
 * @brief       Updates solution array.
 * @param       neles       Number of element cells
 * @param       uptsb       Solution array
 * @param       dub         Result of Block LU-SGS sweeps
 * @param       subres      Residual of each sub-iteration
 */
void blusgs_serial_update(UCFD_INT neles, UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *subres);

#endif // BLUSGS_H