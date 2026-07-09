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
void ns_serial_pre_blusgs(UCFDInt neles, UCFDInt nface, UCFDReal factor,
                          UCFDReal *fnorm_vol, UCFDReal *dt, UCFDReal *diag, UCFDReal *fjmat);


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
void rans_serial_pre_blusgs(UCFDInt neles, UCFDInt nface, UCFDReal factor,
                            UCFDReal *fnorm_vol, UCFDReal *uptsb, UCFDReal *dt,
                            UCFDReal *tdiag, UCFDReal *tjmat, UCFDReal *dsrc);


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
void ns_serial_block_lower_sweep(UCFDInt neles, UCFDInt nface,
                                 UCFDInt *nei_ele, UCFDReal *fnorm_vol,
                                 UCFDReal *rhsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fjmat);


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
void rans_serial_block_lower_sweep(UCFDInt neles, UCFDInt nface,
                                   UCFDInt *nei_ele, UCFDReal *fnorm_vol,
                                   UCFDReal *rhsb, UCFDReal *dub, UCFDReal *tdiag, UCFDReal *tjmat);


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
void ns_serial_block_upper_sweep(UCFDInt neles, UCFDInt nface,
                                 UCFDInt *nei_ele, UCFDReal *fnorm_vol,
                                 UCFDReal *rhsb, UCFDReal *dub, UCFDReal *diag, UCFDReal *fjmat);


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
void rans_serial_block_upper_sweep(UCFDInt neles, UCFDInt nface,
                                   UCFDInt *nei_ele, UCFDReal *fnorm_vol,
                                   UCFDReal *rhsb, UCFDReal *dub, UCFDReal *tdiag, UCFDReal *tjmat);


/**
 * @brief       Updates solution array.
 * @param       neles       Number of element cells
 * @param       uptsb       Solution array
 * @param       dub         Result of Block LU-SGS sweeps
 * @param       subres      Residual of each sub-iteration
 */
void blusgs_serial_ns_update(UCFDInt neles, UCFDReal *uptsb, UCFDReal *dub, UCFDReal *subres);


/**
 * @brief       Updates solution array.
 * @param       neles       Number of element cells
 * @param       uptsb       Solution array
 * @param       dub         Result of Block LU-SGS sweeps
 * @param       subres      Residual of each sub-iteration
 */
void blusgs_serial_update(UCFDInt neles, UCFDReal *uptsb, UCFDReal *dub, UCFDReal *subres);

#endif // BLUSGS_H