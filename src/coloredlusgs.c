/** ======================================================================================================================
 * @file        coloredlusgs.c
 * @brief       UCFD_SPARSE : Unstructured Grid Based CFD Applicable Asymmetric Sparse Matrix Numerical Library
 * @details     Colored LU-SGS time integration method for unstructured grid.  
 *              Multi-thread computation is enabled using <omp.h> header file.
 *              To compute solution of the next time step, refer to the following steps.  
 *                
 *              (1) Preparing LU-SGS :  
 *                  Computes diagonal matrix of the Implicit operator.
 *              (2) Lower sweep :  
 *                  First step of LU-SGS method, which computes intermediate solution difference array.
 *              (3) Upper sweep :  
 *                  Second step of LU-SGS method, which computes next solution difference array.
 *              (4) Update :  
 *                  Time integration by adding current solution with next time step solution difference
 * 
 * @note        In Colored LU-SGS functions, `icolor` and `lcolor` array are used instead of `mapping`, `unmapping` array
 *              which are used in serial LU-SGS. Be aware that they have different values
 *              albeit the number of arguments and each data type is totally identical.
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
 * =======================================================================================================================
 */

#include <omp.h>
#include "coloredlusgs.h"
#include "flux.h"

/**
 * @details     This function computes diagonal matrix of the implicit operator.
 *              In LU-SGS method, implicit operator is replaced with `spectral radius`,
 *              so that all element except diagonal is zero.
 *              Therefore, `diag` array can be allocated as one-dimensional,
 *              which has less memory requirement.  
 *              Diffusive margin of wave speed is applied.
 */
void parallel_pre_lusgs(int neles, int nface, double factor, \
                        double *fnorm_vol, double *dt, double *diag, double *fspr)
{
    int idx;        // Element index
    int jdx;        // Face index
    double lamf;    // Spectral radius at each face

    // SMP applied
    #pragma omp parallel for private(jdx, lamf)
    for (idx=0; idx<neles; idx++) {
        // Diagonals of implicit operator
        diag[idx] = 1.0/(dt[idx]*factor);

        for (jdx=0; jdx<nface; jdx++) {
            // Diffusive margin of wave speed of face
            lamf = fspr[neles*jdx + idx]*1.01;

            // Save spectral radius
            fspr[neles*jdx + idx] = lamf;

            // Add portion of lower and upper spectral radius
            diag[idx] += 0.5*lamf*fnorm_vol[neles*jdx + idx];
        }
    }
}


/**
 * @details     This function computes diagonal matrix of the implicit operator.
 *              In LU-SGS method, implicit operator is replaced with `spectral radius`,
 *              so that all element except diagonal is zero.
 *              Therefore, `diag` array can be allocated as one-dimensional,
 *              which has less memory requirement.  
 *              Diffusive margin of wave speed is applied.
 */
void ns_parallel_lower_sweep(int n0, int ne, int neles, int nfvars, int nface, int ndims, \
                             int *nei_ele, int *icolor, int *lcolor, double *fnorm_vol, double *vec_fnorm, \
                             double *uptsb, double *rhsb, double *dub, double *diag, double *fspr)
{   
    int _idx, idx, jdx, kdx, neib, curr_level;
    double du[nfvars], dfj[nfvars], df[nfvars], nf[ndims];
    double u[nfvars], f[nfvars];
    
    // Lower sweep via coloring
    #pragma omp parallel for private(idx, jdx, kdx, neib, curr_level, \
                                     du, dfj, df, nf, u, f)
    for (_idx=n0; _idx<ne; _idx++) {
        idx = icolor[_idx];
        curr_level = lcolor[idx];

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

            // Only for lower level cell
            if (lcolor[neib] < curr_level) {

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
 * @details     By processing lower sweep, intermediate solution \f$\Delta Q^*$\f is computed.
 *              This function is used for RANS equations.  
 *              solution array is stored in `dub` array.
 */
void rans_parallel_lower_sweep(int n0, int ne, int neles, int nvars, int nfvars, int nface, int ndims, \
                               int *nei_ele, int *icolor, int *lcolor, double *fnorm_vol, double *vec_fnorm, \
                               double *uptsb, double *rhsb, double *dub, double *diag, double *fspr, double *dsrc)
{
    int _idx, idx, jdx, kdx, neib, curr_level;
    int dnv = nvars - nfvars;
    double du[nvars], dfj[dnv], df[dnv], nf[ndims];
    double u[nvars], f[dnv];

    #pragma omp parallel for private(idx, jdx, kdx, neib, curr_level, \
                                     du, dfj, df, nf, u, f)
    // Lower sweep via coloring
    for (_idx=n0; _idx<ne; _idx++) {
        idx = icolor[_idx];
        curr_level = lcolor[idx];

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

            // Only for lower level cell
            if (lcolor[neib] < curr_level) {

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
 * @details     By processing upper sweep, next time step solution \f$\Delta Q$\f is computed.
 *              This function is used for Euler or Navier-Stokes equations,
 *              which has the same flux shape.  
 *              Solution array is stored in `rhsb` array,
 *              since right-hand-side array is no longer needed.
 */
void ns_parallel_upper_sweep(int n0, int ne, int neles, int nfvars, int nface, int ndims, \
                             int *nei_ele, int *icolor, int *lcolor, double *fnorm_vol, double *vec_fnorm, \
                             double *uptsb, double *rhsb, double *dub, double *diag, double *fspr)
{   
    int _idx, idx, jdx, kdx, neib, curr_level;
    double du[nfvars], dfj[nfvars], df[nfvars], nf[ndims];
    double u[nfvars], f[nfvars];

    // Upper sweep via coloring
    #pragma omp parallel for private(idx, jdx, kdx, neib, curr_level, \
                                     du, dfj, df, nf, u, f)
    for (_idx=n0; _idx<ne; _idx++) {
        idx = icolor[_idx];
        curr_level = lcolor[idx];

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

            // Only for upper level cell
            if (lcolor[neib] > curr_level) {

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
 * @details     By processing upper sweep, next time step solution \f$\Delta Q$\f is computed.
 *              This function is used for RANS equations.  
 *              Solution array is stored in `rhsb` array,
 *              since right-hand-side array is no longer needed.
 */
void rans_parallel_upper_sweep(int n0, int ne, int neles, int nvars, int nfvars, int nface, int ndims, \
                               int *nei_ele, int *icolor, int *lcolor, double *fnorm_vol, double *vec_fnorm, \
                               double *uptsb, double *rhsb, double *dub, double *diag, double *fspr, double *dsrc)
{
    int _idx, idx, jdx, kdx, neib, curr_level;
    int dnv = nvars - nfvars;
    double du[nvars], dfj[dnv], df[dnv], nf[ndims];
    double u[nvars], f[dnv];

    #pragma omp parallel for private(idx, jdx, kdx, neib, curr_level, \
                                     du, dfj, df, nf, u, f)
    // Upper sweep via coloring
    for (_idx=n0; _idx<ne; _idx++) {
        idx = icolor[_idx];
        curr_level = lcolor[idx];

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

            // Only for upper level cell
            if (lcolor[neib] > curr_level) {

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
 * @details     solution array updated by adding \f$\Delta Q$\f.
 *              Be aware that `rhsb` array in function parameter
 *              is the difference array after upper sweep,
 *              not the right-hand-side array.
 */
void parallel_update(int neles, int nvars, double *uptsb, double *rhsb)
{
    int idx, kdx;   /** index variables are defined to make private */

    #pragma omp parallel for private(kdx)
    // Iterate for all cell
    for (idx=0; idx<neles; idx++) {
        // Update conservative variables
        for (kdx=0; kdx<nvars; kdx++) {
            // Indexing 2D array as 1D
            uptsb[neles*kdx + idx] += rhsb[neles*kdx + idx];
        }
    }
}
