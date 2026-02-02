/** ======================================================================================================================
 * @file        coloredlusgs.c
 * @brief       Colored LU-SGS time integration method for unstructured grid.
 * @details     LU-SGS time integration method applying Multi-coloring algorithm.  
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
 * 
 * =======================================================================================================================
 */

#include <omp.h>
#include "coloredlusgs.h"

/**
 * @details     This function computes diagonal matrix of the implicit operator.
 *              In LU-SGS method, implicit operator is replaced with `spectral radius`,
 *              so that all element except diagonal is zero.
 *              Therefore, `diag` array can be allocated as one-dimensional,
 *              which has less memory requirement.  
 *              Diffusive margin of wave speed is applied.
 */
void parallel_pre_lusgs(UCFD_INT neles, UCFD_INT nface, UCFD_FLOAT factor,
                        UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *dt, UCFD_FLOAT *diag, UCFD_FLOAT *fspr)
{
    UCFD_INT idx, jdx;
    UCFD_FLOAT lamf;

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
void ns_parallel_lower_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                             UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                             UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr)
{   
    UCFD_INT _idx, idx, jdx, kdx, neib, curr_level;
    UCFD_FLOAT du[NFVARS], dfj[NFVARS], df[NFVARS], nf[NDIMS];
    UCFD_FLOAT u[NFVARS], f[NFVARS];
    
    // Lower sweep via coloring
    #pragma omp parallel for private(idx, jdx, kdx, neib, curr_level, du, dfj, df, nf, u, f)
    for (_idx=n0; _idx<ne; _idx++) {
        idx = icolor[_idx];
        curr_level = lcolor[idx];

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

            // Only for lower level cell
            if (lcolor[neib] < curr_level) {
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
void rans_parallel_lower_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                               UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                               UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr, UCFD_FLOAT *dsrc)
{
    UCFD_INT _idx, idx, jdx, kdx, neib, curr_level;
    UCFD_FLOAT du[NVARS], dfj[NTURBVARS], df[NTURBVARS], nf[NDIMS];
    UCFD_FLOAT u[NVARS], f[NTURBVARS];

    #pragma omp parallel for private(idx, jdx, kdx, neib, curr_level, du, dfj, df, nf, u, f)
    // Lower sweep via coloring
    for (_idx=n0; _idx<ne; _idx++) {
        idx = icolor[_idx];
        curr_level = lcolor[idx];

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

            // Only for lower level cell
            if (lcolor[neib] < curr_level) {
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
void ns_parallel_upper_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                             UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                             UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr)
{   
    UCFD_INT _idx, idx, jdx, kdx, neib, curr_level;
    UCFD_FLOAT du[NFVARS], dfj[NFVARS], df[NFVARS], nf[NDIMS];
    UCFD_FLOAT u[NFVARS], f[NFVARS];

    // Upper sweep via coloring
    #pragma omp parallel for private(idx, jdx, kdx, neib, curr_level, du, dfj, df, nf, u, f)
    for (_idx=n0; _idx<ne; _idx++) {
        idx = icolor[_idx];
        curr_level = lcolor[idx];

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

            // Only for upper level cell
            if (lcolor[neib] > curr_level) {
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
void rans_parallel_upper_sweep(UCFD_INT n0, UCFD_INT ne, UCFD_INT neles, UCFD_INT nface,
                               UCFD_INT *nei_ele, UCFD_INT *icolor, UCFD_INT *lcolor, UCFD_FLOAT *fnorm_vol, UCFD_FLOAT *vec_fnorm,
                               UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb, UCFD_FLOAT *dub, UCFD_FLOAT *diag, UCFD_FLOAT *fspr, UCFD_FLOAT *dsrc)
{
    UCFD_INT _idx, idx, jdx, kdx, neib, curr_level;
    UCFD_FLOAT du[NVARS], dfj[NTURBVARS], df[NTURBVARS], nf[NDIMS];
    UCFD_FLOAT u[NVARS], f[NTURBVARS];

    #pragma omp parallel for private(idx, jdx, kdx, neib, curr_level, du, dfj, df, nf, u, f)
    // Upper sweep via coloring
    for (_idx=n0; _idx<ne; _idx++) {
        idx = icolor[_idx];
        curr_level = lcolor[idx];

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

            // Only for upper level cell
            if (lcolor[neib] > curr_level) {
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
 * @details     solution array is updated by adding \f$\Delta Q\f$.
 *              Be aware that `rhsb` array in function parameter
 *              is the difference array after upper sweep,
 *              not the right-hand-side array.
 */
void clusgs_parallel_update(UCFD_INT neles, UCFD_FLOAT *uptsb, UCFD_FLOAT *rhsb)
{
    UCFD_INT idx, kdx;

    #pragma omp parallel for private(kdx)
    // Iterate for all cell
    for (idx=0; idx<neles; idx++) {
        // Update conservative variables
        for (kdx=0; kdx<NVARS; kdx++) {
            // Indexing 2D array as 1D
            uptsb[neles*kdx + idx] += rhsb[neles*kdx + idx];
        }
    }
}
