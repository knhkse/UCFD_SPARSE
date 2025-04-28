/** ======================================================================================================================
 * @file        coloredblusgs.c
 * @brief       Colored Block LU-SGS time integration method for unstructured grid.
 * @details     Block LU-SGS time integration method with Multi-coloring algorithm.  
 *              Multi-thread computation is enabled using <omp.h> header file.  
 *              In contrast to LU-SGS method, Block LU-SGS method uses `block operator` instead of `spectral radius`.  
 * 
 * @note        In case of RANS equations, add `rans_parallel_pre_blusgs` function right after the `ns_parallel_pre_blusgs` function.  
 *              There is only `block_sweep` function instead of `block_lower_sweep` or `block_upper_sweep` in `blusgs.c`.  
 *              To execute properly, `n0` and `ne` arguments must be input in appropriate order.
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

#include <stdio.h>
#include <omp.h>
#include "coloredblusgs.h"
#include "flux.h"
#include "inverse.h"

/**
 * @details     This function computes diagonal matrices of the implicit operator using multi-thread computing.  
 *              In Block LU-SGS method, implicit operator is approximated with `block operator`.
 *              Diagonal matrices is composed of block operator matrix, which size is n-by-n.
 *              `n` is the number of conservative variables in Navier-Stokes equations.
 */
void ns_parallel_pre_blusgs(int neles, int nfvars, int nface, double factor, \
                            double *fnorm_vol, double *dt, double *diag, double *fjmat)
{
    int idx, jdx, kdx, row, col;        // Element index
    int matsize = nfvars*nfvars;
    double fv, dti;
    double dmat[matsize];     // Diagonal matrix at each cell

    #pragma omp parallel for private(jdx, kdx, row, col, fv, dmat, dti)
    for (idx=0; idx<neles; idx++) {
        // Initialize diagonal matrix
        for (kdx=0; kdx<matsize; kdx++)
            dmat[kdx] = 0.0;
        
        // Computes diagonal matrix based on neighbor cells
        for (jdx=0; jdx<nface; jdx++) {
            fv = fnorm_vol[neles*jdx + idx];
            for (row=0; row<nfvars; row++) {
                for (col=0; col<nfvars; col++) {
                    dmat[nfvars*row+col] \
                        += fjmat[idx+neles*jdx+nface*neles*col+nfvars*nface*neles*row]*fv;
                }
            }
        }

        // Complete implicit operator
        dti = 1.0/(dt[idx]*factor);
        for (kdx=0; kdx<nfvars; kdx++) {
            dmat[(nfvars+1)*kdx] += dti;
        }
        
        // LU decomposition for inverse process
        ludcmp(nfvars, dmat);

        // Allocate temporal matrix to diag array
        for (row=0; row<nfvars; row++) {
            for (col=0; col<nfvars; col++) {
                diag[idx+neles*col+neles*nfvars*row] = dmat[nfvars*row+col];
            }
        }
    }
}


/**
 * @details     This function computes diagonal matrices of the implicit operator of RANS equations.
 *              In Block LU-SGS method, implicit operator is approximated with `block operator`.
 *              Diagonal matrices is composed of block operator matrix, which size is n-by-n.
 *              `n` is the number of turbulent variables in RANS equations.
 */
void rans_parallel_pre_blusgs(int neles, int nvars, int nfvars, int nface, double factor, double betast, \
                              double *fnorm_vol, double *uptsb, double *dt, double *tdiag, double *tjmat, double *dsrc)
{
    int idx, jdx, kdx, row, col;
    int ntvars = nvars - nfvars;        // Constant
    int matsize = ntvars*ntvars;        // Constant
    double tmat[matsize];     // Diagonal matrix at each cell
    double uf[nvars], dsrcf[nvars];
    double fv;
    int err;

    #pragma omp parallel for private(jdx, kdx, row, col, fv, err, \
                                     tmat, uf, dsrcf)
    for (idx=0; idx<neles; idx++) {
        // Initialize diagonal matrix
        for (kdx=0; kdx<matsize; kdx++)
            tmat[kdx] = 0.0;
        
        for (kdx=0; kdx<nvars; kdx++) {
            uf[kdx] = uptsb[idx+neles*kdx];
            dsrcf[kdx] = dsrc[idx+neles*kdx];
        }
        
        // Computes diagonal matrix based on neighbor cells
        for (jdx=0; jdx<nface; jdx++) {
            fv = fnorm_vol[neles*jdx + idx];
            for (row=0; row<ntvars; row++) {
                for (col=0; col<ntvars; col++) {
                    tmat[ntvars*row+col] \
                        += tjmat[idx+neles*jdx+nface*neles*col+ntvars*nface*neles*row]*fv;
                }
            }
        }

        // Computes Source term Jacobian
        err = rans_source_jacobian(nvars, ntvars, betast, uf, tmat, dsrcf);
        if (err == -1) {
            printf("Warning:::Source term Jacobian of RANS equations does not match\n");
        }

        // Complete implicit operator
        for (kdx=0; kdx<ntvars; kdx++) {
            tmat[(ntvars+1)*kdx] += 1.0/(dt[idx]*factor);
        }
        
        // LU decomposition for inverse process
        ludcmp(ntvars, tmat);

        // Allocate temporal matrix to diag array
        for (row=0; row<ntvars; row++) {
            for (col=0; col<ntvars; col++) {
                tdiag[idx+neles*col+neles*ntvars*row] = tmat[ntvars*row+col];
            }
        }
    }
}


/**
 * @details     Lower/Upper sweep of Block LU-SGS with multi-coloring algorithm.  
 *              This function is used for Euler or Navier-Stokes equations,
 *              which has the same flux shape.  
 *              solution array is stored in `dub` array.
 * 
 * @note        The last argument array, `fjmat` is NOT identical with ns_serial_pre_blusgs function.  
 *              For more details, refer to the Block LU-SGS in the document.
 */
void ns_parallel_block_sweep(int n0, int ne, int neles, int nfvars, int nface, \
                             int *nei_ele, int *icolor, int *lcolor, double *fnorm_vol, \
                             double *rhsb, double *dub, double *diag, double *fjmat)
{
    int _idx, idx, jdx, kdx, neib, curr_level;
    int row, col;
    double val, fv;
    double rhs[nfvars], dmat[nfvars*nfvars];

    // Lower/Upper sweep via coloring
    #pragma omp parallel for private(idx, jdx, kdx, neib, curr_level, row, col, \
                                     rhs, dmat, val, fv)
    for (_idx=n0; _idx<ne; _idx++) {
        idx = icolor[_idx];
        curr_level = lcolor[idx];

        // Initialize
        for (kdx=0; kdx<nfvars; kdx++) {
            rhs[kdx] = rhsb[idx+kdx*neles];
        }

        for (row=0; row<nfvars; row++) {
            for (col=0; col<nfvars; col++) {
                dmat[col+nfvars*row] = diag[idx+neles*col+nfvars*neles*row];
            }
        }

        for (jdx=0; jdx<nface; jdx++) {
            neib = nei_ele[idx+neles*jdx];

            if (lcolor[neib] != curr_level) {
                fv = fnorm_vol[idx+neles*jdx];

                for (row=0; row<nfvars; row++) {
                    val = 0.0;
                    for (col=0; col<nfvars; col++) {
                        val += fjmat[idx+neles*jdx+nface*neles*col+nfvars*nface*neles*row] \
                                * dub[neib+neles*col];
                    }
                    rhs[row] -= val*fv;
                }
            }
        }

        lusubst(nfvars, dmat, rhs);

        for (kdx=0; kdx<nfvars; kdx++) {
            dub[idx+neles*kdx] = rhs[kdx];
        }
    }
}


/**
 * @details     Lower/Upper sweep of Block LU-SGS with multi-coloring algorithm.  
 *              This function is used for RANS equations.
 *              solution array is stored in `dub` array.
 * 
 * @note        The last argument array, `tjmat` is NOT identical with `rans_parallel_pre_blusgs` function.  
 *              For more details, refer to the Block LU-SGS in the document.
 */
void rans_parallel_block_sweep(int n0, int ne, int neles, int nvars, int nfvars, int nface, \
                               int *nei_ele, int *icolor, int *lcolor, double *fnorm_vol, \
                               double *rhsb, double *dub, double *tdiag, double *tjmat)
{
    int _idx, idx, jdx, kdx, neib, curr_level;
    int row, col;
    int ntvars = nvars - nfvars;
    double val, fv;
    double rhs[ntvars], dmat[ntvars*ntvars];

    // Lower/Upper sweep via coloring
    #pragma omp parallel for private(idx, jdx, kdx, neib, curr_level, row, col, \
                                     rhs, dmat, val, fv)
    for (_idx=n0; _idx<ne; _idx++) {
        idx = icolor[_idx];
        curr_level = lcolor[idx];

        //Initialize
        for (kdx=0; kdx<ntvars; kdx++) {
            rhs[kdx] = rhsb[idx+(kdx+nfvars)*neles];
        }

        for (row=0; row<ntvars; row++) {
            for (col=0; col<ntvars; col++) {
                dmat[col+ntvars*row] = tdiag[idx+neles*col+neles*ntvars*row];
            }
        }

        for (jdx=0; jdx<nface; jdx++) {
            neib = nei_ele[idx+neles*jdx];

            if (lcolor[neib] != curr_level) {
                fv = fnorm_vol[idx+neles*jdx];
                
                for (row=0; row<ntvars; row++) {
                    val = 0.0;
                    for (col=0; col<ntvars; col++) {
                        val += tjmat[idx+neles*jdx+nface*neles*col+ntvars*nface*neles*row] \
                                * dub[neib+neles*(col+nfvars)];
                    }
                    rhs[row] -= val*fv;
                }
            }
        }

        // Compute inverse of diagonal matrix multiplication
        lusubst(ntvars, dmat, rhs);

        // Update dub array
        for (kdx=0; kdx<ntvars; kdx++) {
            dub[idx+neles*(kdx+nfvars)] = rhs[kdx];
        }
    }
}


/**
 * @details     solution array is updated by adding \f$\Delta Q\f$.
 */
void parallel_update(int neles, int nvars, double *uptsb, double *dub, double *subres)
{
    int idx, kdx;

    #pragma omp parallel for private(kdx)
    for (idx=0; idx<neles; idx++) {
        for (kdx=0; kdx<nvars; kdx++) {
            uptsb[idx+neles*kdx] += dub[idx+neles*kdx];

            // Initialize dub array
            dub[idx+neles*kdx] = 0.0;
        }
        // Initialize sub-residual array
        subres[idx] = 0.0;
    }
}


