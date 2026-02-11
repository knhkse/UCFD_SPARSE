/** ======================================================================================================================
 * @file        blusgs.c
 * @brief       Block LU-SGS time integration method for unstructured grid.
 * @details     Block LU-SGS time integration method (Single thread only).  
 *              In contrast to LU-SGS method, Block LU-SGS method uses `block operator` instead of `spectral radius`.  
 *              Computation process is same with LU-SGS.  
 *              For more information, refer to the `lusgs.c`.  
 * 
 * @note        In case of RANS equations, add rans_serial_pre_blusgs function right after the ns_serial_pre_blusgs function.  
 *              Also, add rans_serial_{}_sweep function right after the ns_serial_{}_sweep function.  
 *              Be aware that ns and rans pre_blusgs/sweep functions must be paired with each sweep step.
 * 
 * @author
 *              - Namhyoung Kim (knhkse@inha.edu), Department of Aerospace Engineering, Inha University
 *              - Jin Seok Park (jinseok.park@inha.ac.kr), Department of Aerospace Engineering, Inha University
 * 
 * @date        Nov 2024
 * @version     1.0
 * @par         Copyright
 *              Copyright (c) 2024, Namhyoung Kim and Jin Seok Park, Inha University, All rights reserved.
 * @par         License
 *              This project is release under the terms of the MIT License (see LICENSE file).
 * 
 * =======================================================================================================================
 */
#include "blusgs.h"
#include "flux.h"
#include "inverse.h"
#include <stdio.h>


/**
 * @details     This function computes diagonal matrices of the implicit operator.
 *              In Block LU-SGS method, implicit operator is approximated with `block operator`.
 *              Diagonal matrices is composed of block operator matrix, which size is n-by-n.
 *              `n` is the number of conservative variables in Navier-Stokes equations,
 *              or the number of turbulent variables in RANS equations.
 */
void ns_serial_pre_blusgs(UCFD_INT neles, UCFD_INT nface, UCFD_FLOAT factor,
                          UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *dt, UCFD_FLOAT *diag, UCFD_FLOAT *fjmat)
{
    UCFD_INT idx, jdx, kdx, row, col;
    UCFD_FLOAT fv, dti;
    UCFD_FLOAT dmat[NFVARS][NFVARS];

    for (idx=0; idx<neles; idx++) {
        // Initialize diagonal matrix
        for (row=0; row<NFVARS; row++) {
            for (col=0; col<NFVARS; col++)
                dmat[row][col] = 0.0;
        }
        
        // Computes diagonal matrix based on neighbor cells
        for (jdx=0; jdx<nface; jdx++) {
            fv = fnorm_vol[neles*jdx + idx];
            for (row=0; row<NFVARS; row++) {
                for (col=0; col<NFVARS; col++) {
                    dmat[row][col] \
                        += fjmat[idx+neles*jdx+nface*neles*col+NFVARS*nface*neles*row]*fv;
                }
            }
        }

        // Complete implicit operator
        dti = 1.0/(dt[idx]*factor);
        for (kdx=0; kdx<NFVARS; kdx++) {
            dmat[kdx][kdx] += dti;
        }
        
        // LU decomposition for inverse process
        ludcmp(NFVARS, dmat);

        // Allocate temporal matrix to diag array
        for (row=0; row<NFVARS; row++) {
            for (col=0; col<NFVARS; col++) {
                diag[idx+neles*col+neles*NFVARS*row] = dmat[row][col];
            }
        }
    }
}


/**
 * @details     This function computes diagonal matrices of the implicit operator of RANS equations.
 *              In Block LU-SGS method, implicit operator is approximated with `block operator`.
 *              Diagonal matrices is composed of block operator matrix, which size is n-by-n.
 *              `n` is the number of conservative variables in Navier-Stokes equations,
 *              or the number of turbulent variables in RANS equations.
 */
void rans_serial_pre_blusgs(UCFD_INT neles, UCFD_INT nface, UCFD_FLOAT factor,
                            UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *uptsb, UCFD_FLOAT *dt,
                            UCFD_FLOAT *tdiag, UCFD_FLOAT *tjmat, UCFD_FLOAT *dsrc)
{
    UCFD_INT idx, jdx, kdx, row, col;
    UCFD_FLOAT fv;
    UCFD_FLOAT tmat[NTURBVARS][NTURBVARS];
    UCFD_FLOAT uf[NVARS], dsrcf[NVARS];

    for (idx=0; idx<neles; idx++) {
        // Initialize diagonal matrix
        for (row=0; row<NTURBVARS; row++) {
            for (col=0; col<NTURBVARS; col++)
                tmat[row][col] = 0.0;
        }
        
        for (kdx=0; kdx<NVARS; kdx++) {
            uf[kdx] = uptsb[idx+neles*kdx];
            dsrcf[kdx] = dsrc[idx+neles*kdx];
        }
        
        // Computes diagonal matrix based on neighbor cells
        for (jdx=0; jdx<nface; jdx++) {
            fv = fnorm_vol[neles*jdx + idx];
            for (row=0; row<NTURBVARS; row++) {
                for (col=0; col<NTURBVARS; col++) {
                    tmat[row][col] \
                        += tjmat[idx+neles*jdx+nface*neles*col+NTURBVARS*nface*neles*row]*fv;
                }
            }
        }

        // Computes Source term Jacobian
        if (rans_source_jacobian(uf, tmat, dsrcf) == UCFD_STATUS_NOT_SUPPORTED)
            printf("Error::Invalid `NTURBVARS` value\n");

        // Complete implicit operator
        for (kdx=0; kdx<NTURBVARS; kdx++) {
            tmat[kdx][kdx] += 1.0/(dt[idx]*factor);
        }
        
        // LU decomposition for inverse process
        ludcmp(NTURBVARS, tmat);

        // Allocate temporal matrix to diag array
        for (row=0; row<NTURBVARS; row++) {
            for (col=0; col<NTURBVARS; col++) {
                tdiag[idx+neles*col+neles*NTURBVARS*row] = tmat[row][col];
            }
        }
    }
}


/**
 * @details     By processing lower sweep, intermediate solution \f$\Delta Q^*\f$ is computed.
 *              This function is used for Euler or Navier-Stokes equations,
 *              which has the same flux shape.  
 *              solution array is stored in `dub` array.
 * 
 * @note        The last argument array, `fjmat` is NOT identical with ns_serial_pre_blusgs function.  
 *              For more details, refer to the Block LU-SGS in the document.
 */
void ns_serial_block_lower_sweep(UCFD_INT neles, UCFD_INT nface,
                                 UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol,
                                 UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fjmat)
{
    UCFD_INT idx, jdx, kdx, neib;
    UCFD_INT row, col;
    UCFD_FLOAT rhs[NFVARS], dmat[NFVARS][NFVARS];
    UCFD_FLOAT val, fv;

    // Lower sweep via mapping
    for (idx=0; idx<neles; idx++) {
        // Initialize
        for (kdx=0; kdx<NFVARS; kdx++) {
            rhs[kdx] = rhsb[idx+kdx*neles];
        }

        for (row=0; row<NFVARS; row++) {
            for (col=0; col<NFVARS; col++) {
                dmat[row][col] = diag[idx+neles*col+neles*NFVARS*row];
            }
        }

        // Only for neighbor cells
        for (jdx=0; jdx<nface; jdx++) {
            neib = nei_ele[idx+neles*jdx];

            if (neib != idx) {
                fv = fnorm_vol[neles*jdx + idx];
                // Matrix-Vector multiplication
                for (row=0; row<NFVARS; row++) {
                    val = 0.0;
                    for (col=0; col<NFVARS; col++) {
                        val += fjmat[idx+neles*jdx+nface*neles*col+NFVARS*nface*neles*row] \
                                * dub[neib+neles*col];
                    }
                    rhs[row] -= val*fv;
                }
            }
        }

        // Compute inverse of diagonal matrix multiplication
        lusub(NFVARS, dmat, rhs);

        // Update dub array
        for (kdx=0; kdx<NFVARS; kdx++) {
            dub[idx+neles*kdx] = rhs[kdx];
        }
    }
}


/**
 * @details     Lower sweep of Block LU-SGS.  
 *              This function is used for RANS equations.
 *              solution array is stored in `dub` array.
 * 
 * @note        The last argument array, `tjmat` is NOT identical with ns_serial_pre_blusgs function.  
 *              For more details, refer to the Block LU-SGS in the document.
 */
void rans_serial_block_lower_sweep(UCFD_INT neles, UCFD_INT nface,
                                   UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol,
                                   UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *tdiag, UCFD_FLOAT *tjmat)
{
    UCFD_INT idx, jdx, kdx, neib;
    UCFD_INT row, col;
    UCFD_FLOAT val, fv;
    UCFD_FLOAT rhs[NTURBVARS], tmat[NTURBVARS][NTURBVARS];

    // Lower sweep via mapping
    for (idx=0; idx<neles; idx++) {
        // Initialize
        for (kdx=0; kdx<NTURBVARS; kdx++) {
            rhs[kdx] = rhsb[idx+(kdx+NFVARS)*neles];
        }

        for (row=0; row<NTURBVARS; row++) {
            for (col=0; col<NTURBVARS; col++) {
                tmat[row][col] = tdiag[idx+neles*col+neles*NTURBVARS*row];
            }
        }

        // Only for neighbor cells
        for (jdx=0; jdx<nface; jdx++) {
            neib = nei_ele[idx+neles*jdx];

            if (neib != idx) {
                fv = fnorm_vol[idx+neles*jdx];
                // Matrix-Vector multiplication
                for (row=0; row<NTURBVARS; row++) {
                    val = 0.0;
                    for (col=0; col<NTURBVARS; col++) {
                        val += tjmat[idx+neles*jdx+nface*neles*col+NTURBVARS*nface*neles*row] \
                                * dub[neib+neles*(col+NFVARS)];
                    }
                    rhs[row] -= val*fv;
                }
            }
        }

        // Compute inverse of diagonal matrix multiplication
        lusub(NTURBVARS, tmat, rhs);

        // Update dub array
        for (kdx=0; kdx<NTURBVARS; kdx++) {
            dub[idx+neles*(kdx+NFVARS)] = rhs[kdx];
        }
    }
}


/**
 * @details     By processing upper sweep, next sub-iteration solution \f$\Delta Q^(k+1)\f$ is computed.
 *              This function is used for Euler or Navier-Stokes equations,
 *              which has the same flux shape.  
 *              solution array is stored in `dub` array.
 * 
 * @note        The last argument array, `fjmat` is NOT identical with ns_serial_pre_blusgs function.  
 *              For more details, refer to the Block LU-SGS in the document.
 */
void ns_serial_block_upper_sweep(UCFD_INT neles, UCFD_INT nface,
                                 UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol,
                                 UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fjmat)
{
    UCFD_INT idx, jdx, kdx, neib;
    UCFD_INT row, col;
    UCFD_FLOAT val, fv;
    UCFD_FLOAT rhs[NFVARS], dmat[NFVARS][NFVARS];

    // Upper sweep via mapping
    for (idx=neles-1; idx>-1; idx--) {
        // Initialize
        for (kdx=0; kdx<NFVARS; kdx++) {
            rhs[kdx] = rhsb[idx+kdx*neles];
        }

        for (row=0; row<NFVARS; row++) {
            for (col=0; col<NFVARS; col++) {
                dmat[row][col] = diag[idx+neles*col+neles*NFVARS*row];
            }
        }

        // Only for neighbor cells
        for (jdx=0; jdx<nface; jdx++) {
            neib = nei_ele[idx+neles*jdx];

            if (neib != idx) {
                fv = fnorm_vol[neles*jdx + idx];
                // Matrix-Vector multiplication
                for (row=0; row<NFVARS; row++) {
                    val = 0.0;
                    for (col=0; col<NFVARS; col++) {
                        val += fjmat[idx+neles*jdx+nface*neles*col+NFVARS*nface*neles*row] \
                                * dub[neib+neles*col];
                    }
                    rhs[row] -= val*fv;
                }
            }
        }

        // Compute inverse of diagonal matrix multiplication
        lusub(NFVARS, dmat, rhs);

        // Update dub array
        for (kdx=0; kdx<NFVARS; kdx++) {
            dub[idx+neles*kdx] = rhs[kdx];
        }
    }
}


/**
 * @details     Upper sweep of Block LU-SGS.  
 *              This function is used for RANS equations.
 *              solution array is stored in `dub` array.
 * 
 * @note        The last argument array, `tjmat` is NOT identical with ns_serial_pre_blusgs function.  
 *              For more details, refer to the Block LU-SGS in the document.
 */
void rans_serial_block_upper_sweep(UCFD_INT neles, UCFD_INT nface,
                                   UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol,
                                   UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *tdiag, UCFD_FLOAT *tjmat)
{
    UCFD_INT idx, jdx, kdx, neib;
    UCFD_INT row, col;
    UCFD_FLOAT val, fv;
    UCFD_FLOAT rhs[NTURBVARS], tmat[NTURBVARS][NTURBVARS];

    // Lower sweep via mapping
    for (idx=neles-1; idx>-1; idx--) {

        // Initialize
        for (kdx=0; kdx<NTURBVARS; kdx++) {
            rhs[kdx] = rhsb[idx+(kdx+NFVARS)*neles];
        }

        for (row=0; row<NTURBVARS; row++) {
            for (col=0; col<NTURBVARS; col++) {
                tmat[row][col] = tdiag[idx+neles*col+neles*NTURBVARS*row];
            }
        }

        // Only for neighbor cells
        for (jdx=0; jdx<nface; jdx++) {
            neib = nei_ele[idx+neles*jdx];

            if (neib != idx) {
                fv = fnorm_vol[neles*jdx + idx];
                // Matrix-Vector multiplication
                for (row=0; row<NTURBVARS; row++) {
                    val = 0.0;
                    for (col=0; col<NTURBVARS; col++) {
                        val += tjmat[idx+neles*jdx+nface*neles*col+NTURBVARS*nface*neles*row] \
                                * dub[neib+neles*(col+NFVARS)];
                    }
                    rhs[row] -= val*fv;
                }
            }
        }

        // Compute inverse of diagonal matrix multiplication
        lusub(NTURBVARS, tmat, rhs);

        // Update dub array
        for (kdx=0; kdx<NTURBVARS; kdx++) {
            dub[idx+neles*(kdx+NFVARS)] = rhs[kdx];
        }
    }
}


/**
 * @details     solution array is updated by adding \f$\Delta Q\f$.
 */
void blusgs_serial_ns_update(UCFD_INT neles, UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *subres)
{
    UCFD_INT idx, kdx;

    for (idx=0; idx<neles; idx++) {
        for (kdx=0; kdx<NFVARS; kdx++) {
            uptsb[idx+neles*kdx] += dub[idx+neles*kdx];

            // Initialize dub array
            dub[idx+neles*kdx] = 0.0;
        }
        // Initialize sub-residual array
        subres[idx] = 0.0;
    }
}


/**
 * @details     solution array is updated by adding \f$\Delta Q\f$.
 */
void blusgs_serial_update(UCFD_INT neles, UCFD_FLOAT *uptsb, UCFD_FLOAT *dub, UCFD_FLOAT *subres)
{
    UCFD_INT idx, kdx;

    for (idx=0; idx<neles; idx++) {
        for (kdx=0; kdx<NVARS; kdx++) {
            uptsb[idx+neles*kdx] += dub[idx+neles*kdx];

            // Initialize dub array
            dub[idx+neles*kdx] = 0.0;
        }
        // Initialize sub-residual array
        subres[idx] = 0.0;
    }
}