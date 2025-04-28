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
#include "flux.h"

/**
 * @details     This function computes diagonal matrix of the implicit operator.
 *              In LU-SGS method, implicit operator is replaced with `spectral radius`,
 *              so that all element except diagonal is zero.
 *              Therefore, `diag` array can be allocated as one-dimensional,
 *              which has less memory requirement.  
 *              Diffusive margin of wave speed is applied.
 */
void serial_pre_lusgs(int neles, int nface, double factor, \
                      double *fnorm_vol, double *dt, double *diag, double *fspr)
{
    int idx;        // Element index
    int jdx;        // Face index
    double lamf;    // Spectral radius at each face

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
void ns_serial_lower_sweep(int neles, int nfvars, int nface, int ndims, \
                           int *nei_ele, int *mapping, int *unmapping, double *fnorm_vol, double *vec_fnorm, \
                           double *uptsb, double *rhsb, double *dub, double *diag, double *fspr)
{   
    int _idx, idx, jdx, kdx, neib;                              // Index variables
    double du[nfvars], dfj[nfvars], df[nfvars], nf[ndims];      // Temporal arrays
    double u[nfvars], f[nfvars];
    
    // Lower sweep via mapping
    for (_idx=0; _idx<neles; _idx++) {
        idx = mapping[_idx];

        // Initialize `df` array
        for (kdx=0; kdx<nfvars; kdx++) {
            df[kdx] = 0.0;
        }
        
        // Set of faces surrounding a cell
        for (jdx=0; jdx<nface; jdx++) {
            
            // Get face normal vector
            for (kdx=0; kdx<ndims; kdx++) {
                nf[kdx] = vec_fnorm[ndims*neles*jdx + neles*kdx + idx];
            }

            // Neighbor element index meet at face
            neib = nei_ele[neles*jdx + idx];

            // Only for lower neighbor cell
            if (unmapping[neib] < _idx) {

                for (kdx=0; kdx<nfvars; kdx++) {
                    u[kdx] = uptsb[neles*kdx + neib];
                    du[kdx] = u[kdx] + dub[neles*kdx + neib];
                }
                
                ns_flux_container(nfvars, ndims, u, nf, f);
                ns_flux_container(nfvars, ndims, du, nf, dfj);

                for (kdx=0; kdx<nfvars; kdx++) {
                    dfj[kdx] -= f[kdx];
                }

                for (kdx=0; kdx<nfvars; kdx++) {
                    df[kdx] += (dfj[kdx] - fspr[neles*jdx + idx] \
                                * dub[neles*kdx + neib])*fnorm_vol[neles*jdx + idx];
                }
            }
        }

        // Update dub array
        for (kdx=0; kdx<nfvars; kdx++)
            dub[neles*kdx + idx] = (rhsb[neles*kdx + idx] - 0.5*df[kdx])/diag[idx];
    }
}

/**
 * @details     By processing lower sweep, intermediate solution \f$\Delta Q^*\f$ is computed.
 *              This function is used for RANS equations.  
 *              solution array is stored in `dub` array.
 */
void rans_serial_lower_sweep(int neles, int nvars, int nfvars, int nface, int ndims, \
                             int *nei_ele, int *mapping, int *unmapping, double *fnorm_vol, double *vec_fnorm, \
                             double *uptsb, double *rhsb, double *dub, double *diag, double *fspr, double *dsrc)
{
    int _idx, idx, jdx, kdx, neib;
    int dnv = nvars - nfvars;
    double du[nvars], dfj[dnv], df[dnv], nf[ndims];
    double u[nvars], f[dnv];

    // Lower sweep via mapping
    for (_idx=0; _idx<neles; _idx++) {
        idx = mapping[_idx];

        // Initialize `df` array
        for (kdx=0; kdx<dnv; kdx++) {
            df[kdx] = 0.0;
        }
        
        // Set of faces surrounding a cell
        for (jdx=0; jdx<nface; jdx++) {

            // Get face normal vector
            for (kdx=0; kdx<ndims; kdx++) {
                nf[kdx] = vec_fnorm[ndims*neles*jdx + neles*kdx + idx];
            }
            
            // Neighbor element index meet at face
            neib = nei_ele[neles*jdx + idx];

            // Only for lower neighbor cell
            if (unmapping[neib] < _idx) {

                for (kdx=0; kdx<nvars; kdx++) {
                    u[kdx] = uptsb[neles*kdx + neib];
                    du[kdx] = u[kdx];
                }

                for (kdx=nfvars; kdx<nvars; kdx++) {
                    du[kdx] += dub[neles*kdx + neib];
                }
                
                rans_flux_container(nfvars, ndims, dnv, u, nf, f);
                rans_flux_container(nfvars, ndims, dnv, du, nf, dfj);

                for (kdx=0; kdx<dnv; kdx++) {
                    dfj[kdx] -= f[kdx];
                }

                for (kdx=0; kdx<dnv; kdx++) {
                    df[kdx] += (dfj[kdx] - fspr[neles*jdx + idx] \
                                * dub[neles*(kdx+nfvars) + neib]) * fnorm_vol[neles*jdx + idx];
                }
            }
        }

        // Update dub array
        for (kdx=0; kdx<dnv; kdx++) {
            dub[neles*(kdx+nfvars) + idx] = (rhsb[neles*(kdx+nfvars)+idx] - \
                                            0.5*df[kdx])/(diag[idx]+dsrc[neles*(kdx+nfvars)+idx]);
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
void ns_serial_upper_sweep(int neles, int nfvars, int nface, int ndims, \
                           int *nei_ele, int *mapping, int *unmapping, double *fnorm_vol, double *vec_fnorm, \
                           double *uptsb, double *rhsb, double *dub, double *diag, double *fspr)
{
    int _idx, idx, jdx, kdx, neib;
    double du[nfvars], dfj[nfvars], df[nfvars], nf[ndims];
    double u[nfvars], f[nfvars];
    
    // Upper sweep via mapping
    for (_idx=neles-1; _idx>-1; _idx--) {
        idx = mapping[_idx];

        // Initialize `df` array
        for (kdx=0; kdx<nfvars; kdx++) {
            df[kdx] = 0.0;
        }
        
        // Set of faces surrounding a cell
        for (jdx=0; jdx<nface; jdx++) {
            
            // Get face normal vector
            for (kdx=0; kdx<ndims; kdx++) {
                nf[kdx] = vec_fnorm[ndims*neles*jdx + neles*kdx + idx];
            }
            
            // Neighbor element index meet at face
            neib = nei_ele[neles*jdx + idx];

            // Only for upper neighbor cell
            if (unmapping[neib] > _idx) {

                for (kdx=0; kdx<nfvars; kdx++) {
                    u[kdx] = uptsb[neles*kdx + neib];
                    du[kdx] = u[kdx] + rhsb[neles*kdx + neib];
                }
                
                ns_flux_container(nfvars, ndims, u, nf, f);
                ns_flux_container(nfvars, ndims, du, nf, dfj);

                for (kdx=0; kdx<nfvars; kdx++) {
                    dfj[kdx] -= f[kdx];
                }
                
                for (kdx=0; kdx<nfvars; kdx++) {
                    df[kdx] += (dfj[kdx] - fspr[neles*jdx + idx] \
                                * rhsb[neles*kdx + neib])*fnorm_vol[neles*jdx + idx];
                }
            }
        }

        // Update rhsb array
        for (kdx=0; kdx<nfvars; kdx++) {
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
void rans_serial_upper_sweep(int neles, int nvars, int nfvars, int nface, int ndims, \
                             int *nei_ele, int *mapping, int *unmapping, double *fnorm_vol, double *vec_fnorm, \
                             double *uptsb, double *rhsb, double *dub, double *diag, double *fspr, double *dsrc)
{
    int _idx, idx, jdx, kdx, neib;
    int dnv = nvars - nfvars;
    double du[nvars], dfj[dnv], df[dnv], nf[ndims];
    double u[nvars], f[dnv];

    // Upper sweep via mapping
    for (_idx=neles-1; _idx>-1; _idx--) {
        idx = mapping[_idx];

        // Initialize `df` array
        for (kdx=0; kdx<dnv; kdx++) {
            df[kdx] = 0.0;
        }
        
        // Set of faces surrounding a cell
        for (jdx=0; jdx<nface; jdx++) {

            // Get face normal vector
            for (kdx=0; kdx<ndims; kdx++) {
                nf[kdx] = vec_fnorm[ndims*neles*jdx + neles*kdx + idx];
            }

            // Neighbor element index meet at face
            neib = nei_ele[neles*jdx + idx];

            // Only for upper neighbor cell
            if (unmapping[neib] > _idx) {

                for (kdx=0; kdx<nvars; kdx++) {
                    u[kdx] = uptsb[neles*kdx + neib];
                    du[kdx] = u[kdx];
                }

                for (kdx=nfvars; kdx<nvars; kdx++) {
                    du[kdx] += rhsb[neles*kdx + neib];
                }

                rans_flux_container(nfvars, ndims, dnv, u, nf, f);
                rans_flux_container(nfvars, ndims, dnv, du, nf, dfj);

                for (kdx=0; kdx<dnv; kdx++) {
                    dfj[kdx] -= f[kdx];
                }

                for (kdx=0; kdx<dnv; kdx++) {
                    df[kdx] += (dfj[kdx] - fspr[neles*jdx+idx] \
                                * rhsb[neles*(kdx+nfvars)+neib])*fnorm_vol[neles*jdx+idx];
                }
            }
        }

        // Update rhsb array
        for (kdx=0; kdx<dnv; kdx++) {
            rhsb[neles*(kdx+nfvars)+idx] = dub[neles*(kdx+nfvars)+idx] - \
                                        0.5*df[kdx]/(diag[idx] + dsrc[neles*(kdx+nfvars)+idx]);
        }
    }
}


/**
 * @details     solution array updated by adding \f$\Delta Q\f$.
 *              Be aware that `rhsb` array in function parameter
 *              is the difference array after upper sweep,
 *              not the right-hand-side array.
 */
void serial_update(int neles, int nvars, double *uptsb, double *rhsb)
{
    // Iterate for all cell
    for (int idx=0; idx<neles; idx++) {
        // Update conservative variables
        for (int kdx=0; kdx<nvars; kdx++) {
            // Indexing 2D array as 1D
            uptsb[neles*kdx + idx] += rhsb[neles*kdx + idx];
        }
    }
}
