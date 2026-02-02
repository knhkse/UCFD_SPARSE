/** ======================================================================================================================
 * @file        lusgs.c
 * @brief       LU-SGS time integration method for unstructured grid.
 * @details     LU-SGS time integration method (Single thread only).  
 *              To compute solution of the next time step, refer to the following steps.  
 *                
 *              (1) Preparing LU-SGS :  
 *                  Computes diagonal matrix of the Implicit operator.  
 *              (2) Lower sweep :  
 *                  First step of LU-SGS method, which computes intermediate solution difference array.  
 *              (3) Upper sweep :  
 *                  Second step of LU-SGS method, which computes next solution difference array.  
 *              (4) Update :  
 *                  Time integration by adding current solution with next time step solution difference.  
 * 
 * @note        In case of RANS equations, add rans_serial_{}_sweep function right after the ns_serial_{}_sweep function.  
 *              Be aware that ns and rans sweep function must be paired with each sweep step.
 * 
 * @author
 *              - Namhyoung Kim (knhkse@inha.edu), Department of Aerospace Engineering, Inha University
 *              - Jin Seok Park (jinseok.park@inha.ac.kr), Department of Aerospace Engineering, Inha University
 * 
 * @date        July 2024
 * @version     1.0
 * @par         Copyright
 *              Copyright (c) 2024, Namhyoung Kim and Jin Seok Park, Inha University, All rights reserved.
 * @par         License
 *              This project is release under the terms of the MIT License (see LICENSE file).
 * 
 * =======================================================================================================================
 */

#include "lusgs.h"

/**
 * @details     This function computes diagonal matrix of the implicit operator.
 *              In LU-SGS method, implicit operator is replaced with `spectral radius`,
 *              so that all element except diagonal is zero.
 *              Therefore, `diag` array can be allocated as one-dimensional,
 *              which has less memory requirement.  
 *              Diffusive margin of wave speed is applied.
 */
void serial_pre_lusgs(UCFD_INT neles, UCFD_INT nface, UCFD_FLOAT factor, \
                      UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *dt, UCFD_FLOAT *diag, UCFD_FLOAT *fspr)
{
    UCFD_INT idx;       // Element index
    UCFD_INT jdx;       // Face index
    UCFD_FLOAT lamf;    // Spectral radius at each face

    for (idx=0; idx<neles; idx++) {
        // Diagonals of implicit operator
        diag[idx] = 1.0/(dt[idx]*factor);

        for (jdx=0; jdx<nface; jdx++) {
            // Apply diffusive margin of wave speed at face
            lamf = fspr[neles*jdx + idx]*1.01;

            // Save spectral radius
            fspr[neles*jdx + idx] = lamf;

            // Add portion of lower and upper spectral radius
            diag[idx] += 0.5*lamf*fnorm_vol[neles*jdx + idx];
        }
    }
}


/**
 * @details     By processing lower sweep, intermediate solution \f$\Delta Q^*\f$ is computed.
 *              This function is used for Euler or Navier-Stokes equations,
 *              which has the same flux shape.  
 *              solution array is stored in `dub` array.
 */
void ns_serial_lower_sweep(UCFD_INT neles, UCFD_INT nface, UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                           UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr)
{   
    UCFD_INT idx, neib, jdx, kdx;
    UCFD_FLOAT du[NFVARS], dfj[NFVARS], df[NFVARS], nf[NDIMS];
    UCFD_FLOAT u[NFVARS], f[NFVARS];
    
    // Lower sweep via mapping
    for (idx=0; idx<neles; idx++) {
        // Initialize `df` array
        for (kdx=0; kdx<NFVARS; kdx++) {
            df[kdx] = 0.0;
        }

        // Set of faces surrounding a cell
        for (jdx=0; jdx<nface; jdx++) {
            // Get face normal vector
            for (kdx=0; kdx<NDIMS; kdx++) {
                nf[kdx] = vec_fnorm[NDIMS*neles*jdx + neles*kdx + idx];
            }

            // Neighbor element index meet at face
            neib = nei_ele[neles*jdx + idx];

            // Only for lower neighbor cell
            if (neib < idx) {
                for (kdx=0; kdx<NFVARS; kdx++) {
                    u[kdx] = uptsb[neles*kdx + neib];
                    du[kdx] = u[kdx] + dub[neles*kdx + neib];
                }
                
                ns_flux_container(u, nf, f);
                ns_flux_container(du, nf, dfj);

                for (kdx=0; kdx<NFVARS; kdx++) {
                    dfj[kdx] -= f[kdx];
                }

                for (kdx=0; kdx<NFVARS; kdx++) {
                    df[kdx] += (dfj[kdx] - fspr[neles*jdx + idx] \
                                * dub[neles*kdx + neib])*fnorm_vol[neles*jdx + idx];
                }
            }
        }
        // Update dub array
        for (kdx=0; kdx<NFVARS; kdx++)
            dub[neles*kdx + idx] = (rhsb[neles*kdx + idx] - 0.5*df[kdx])/diag[idx];
    }
}

/**
 * @details     By processing lower sweep, intermediate solution \f$\Delta Q^*\f$ is computed.
 *              This function is used for RANS equations.  
 *              solution array is stored in `dub` array.
 */
void rans_serial_lower_sweep(UCFD_INT neles, UCFD_INT nface, UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm, 
                             UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr, UCFD_FLOAT *dsrc)
{
    UCFD_INT idx, neib, jdx, kdx;
    UCFD_FLOAT du[NVARS], dfj[NTURBVARS], df[NTURBVARS], nf[NDIMS];
    UCFD_FLOAT u[NVARS], f[NTURBVARS];

    // Lower sweep via mapping
    for (idx=0; idx<neles; idx++) {
        // Initialize `df` array
        for (kdx=0; kdx<NTURBVARS; kdx++) {
            df[kdx] = 0.0;
        }
        
        // Set of faces surrounding a cell
        for (jdx=0; jdx<nface; jdx++) {
            // Get face normal vector
            for (kdx=0; kdx<NDIMS; kdx++) {
                nf[kdx] = vec_fnorm[NDIMS*neles*jdx + neles*kdx + idx];
            }
            
            // Neighbor element index meet at face
            neib = nei_ele[neles*jdx + idx];

            // Only for lower neighbor cell
            if (neib < idx) {
                for (kdx=0; kdx<NVARS; kdx++) {
                    u[kdx] = uptsb[neles*kdx + neib];
                    du[kdx] = u[kdx];
                }

                for (kdx=NFVARS; kdx<NVARS; kdx++) {
                    du[kdx] += dub[neles*kdx + neib];
                }
                
                rans_flux_container(u, nf, f);
                rans_flux_container(du, nf, dfj);

                for (kdx=0; kdx<NTURBVARS; kdx++) {
                    dfj[kdx] -= f[kdx];
                }

                for (kdx=0; kdx<NTURBVARS; kdx++) {
                    df[kdx] += (dfj[kdx] - fspr[neles*jdx + idx] \
                                * dub[neles*(kdx+NFVARS) + neib]) * fnorm_vol[neles*jdx + idx];
                }
            }
        }
        // Update dub array
        for (kdx=0; kdx<NTURBVARS; kdx++) {
            dub[neles*(kdx+NFVARS) + idx] = (rhsb[neles*(kdx+NFVARS)+idx] - \
                                            0.5*df[kdx])/(diag[idx]+dsrc[neles*(kdx+NFVARS)+idx]);
        }
    }
}


/**
 * @details     By processing upper sweep, next time step solution \f$\Delta Q\f$ is computed.
 *              This function is used for Euler or Navier-Stokes equations,
 *              which has the same flux shape.  
 *              Solution array is stored in `rhsb` array,
 *              since right-hand-side array is no longer needed.
 */
void ns_serial_upper_sweep(UCFD_INT neles, UCFD_INT nface, UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                           UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr)
{
    UCFD_INT idx, neib, jdx, kdx;
    UCFD_FLOAT du[NFVARS], dfj[NFVARS], df[NFVARS], nf[NDIMS];
    UCFD_FLOAT u[NFVARS], f[NFVARS];
    
    // Upper sweep via mapping
    for (idx=neles-1; idx>-1; idx--) {
        // Initialize `df` array
        for (kdx=0; kdx<NFVARS; kdx++) {
            df[kdx] = 0.0;
        }
        
        // Set of faces surrounding a cell
        for (jdx=0; jdx<nface; jdx++) {
            // Get face normal vector
            for (kdx=0; kdx<NDIMS; kdx++) {
                nf[kdx] = vec_fnorm[NDIMS*neles*jdx + neles*kdx + idx];
            }
            
            // Neighbor element index meet at face
            neib = nei_ele[neles*jdx + idx];

            // Only for upper neighbor cell
            if (neib > idx) {
                for (kdx=0; kdx<NFVARS; kdx++) {
                    u[kdx] = uptsb[neles*kdx + neib];
                    du[kdx] = u[kdx] + rhsb[neles*kdx + neib];
                }
                
                ns_flux_container(u, nf, f);
                ns_flux_container(du, nf, dfj);

                for (kdx=0; kdx<NFVARS; kdx++) {
                    dfj[kdx] -= f[kdx];
                }
                
                for (kdx=0; kdx<NFVARS; kdx++) {
                    df[kdx] += (dfj[kdx] - fspr[neles*jdx + idx] \
                                * rhsb[neles*kdx + neib])*fnorm_vol[neles*jdx + idx];
                }
            }
        }
        // Update rhsb array
        for (kdx=0; kdx<NFVARS; kdx++) {
            rhsb[neles*kdx + idx] = dub[neles*kdx + idx] - 0.5*df[kdx]/diag[idx];
        }
    }
}


/**
 * @details     By processing upper sweep, next time step solution \f$\Delta Q\f$ is computed.
 *              This function is used for RANS equations.  
 *              Solution array is stored in `rhsb` array,
 *              since right-hand-side array is no longer needed.
 */
void rans_serial_upper_sweep(UCFD_INT neles, UCFD_INT nface, UCFD_INT *nei_ele, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                             UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr, UCFD_FLOAT *dsrc)
{
    UCFD_INT idx, neib, jdx, kdx;
    UCFD_FLOAT du[NVARS], dfj[NTURBVARS], df[NTURBVARS], nf[NDIMS];
    UCFD_FLOAT u[NVARS], f[NTURBVARS];

    // Upper sweep via mapping
    for (idx=neles-1; idx>-1; idx--) {
        // Initialize `df` array
        for (kdx=0; kdx<NTURBVARS; kdx++) {
            df[kdx] = 0.0;
        }
        
        // Set of faces surrounding a cell
        for (jdx=0; jdx<nface; jdx++) {
            // Get face normal vector
            for (kdx=0; kdx<NDIMS; kdx++) {
                nf[kdx] = vec_fnorm[NDIMS*neles*jdx + neles*kdx + idx];
            }

            // Neighbor element index meet at face
            neib = nei_ele[neles*jdx + idx];

            // Only for upper neighbor cell
            if (neib > idx) {
                for (kdx=0; kdx<NVARS; kdx++) {
                    u[kdx] = uptsb[neles*kdx + neib];
                    du[kdx] = u[kdx];
                }

                for (kdx=NFVARS; kdx<NVARS; kdx++) {
                    du[kdx] += rhsb[neles*kdx + neib];
                }

                rans_flux_container(u, nf, f);
                rans_flux_container(du, nf, dfj);

                for (kdx=0; kdx<NTURBVARS; kdx++) {
                    dfj[kdx] -= f[kdx];
                }

                for (kdx=0; kdx<NTURBVARS; kdx++) {
                    df[kdx] += (dfj[kdx] - fspr[neles*jdx+idx] \
                                * rhsb[neles*(kdx+NFVARS)+neib])*fnorm_vol[neles*jdx+idx];
                }
            }
        }
        // Update rhsb array
        for (kdx=0; kdx<NTURBVARS; kdx++) {
            rhsb[neles*(kdx+NFVARS)+idx] = dub[neles*(kdx+NFVARS)+idx] - \
                                        0.5*df[kdx]/(diag[idx] + dsrc[neles*(kdx+NFVARS)+idx]);
        }
    }
}


/**
 * @details     solution array updated by adding \f$\Delta Q\f$.
 *              Be aware that `rhsb` array in function parameter
 *              is the difference array after upper sweep,
 *              not the right-hand-side array.
 */
void lusgs_serial_update(UCFD_INT neles, UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb)
{
    UCFD_INT idx, kdx;
    
    // Iterate for all cell
    for (idx=0; idx<neles; idx++) {
        // Update conservative variables
        for (kdx=0; kdx<NVARS; kdx++) {
            // Indexing 2D array as 1D
            uptsb[neles*kdx + idx] += rhsb[neles*kdx + idx];
        }
    }
}
