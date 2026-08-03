/** ======================================================================================================================
 * @file        flux.cuh
 * @brief       Computes numerical flux from conservative variables.
 * @details     Non-physical value correction included.
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
#ifndef FLUX_CUDA_H
#define FLUX_CUDA_H
#include "ucfdtypes.h"
#include "config.h"


// Single precision
#if defined(UCFD_FLOAT32)
    #ifndef BETAST
        #define BETAST 0.09f
    #endif
    #ifndef GAMMA
        #define GAMMA 1.4f
    #endif
    #ifndef PMIN
        #define PMIN 1e-13f
    #endif
// Double precision
#else
    #ifndef BETAST
        #define BETAST 0.09
    #endif
    #ifndef GAMMA
        #define GAMMA 1.4
    #endif
    #ifndef PMIN
        #define PMIN 1e-13
    #endif
#endif

/**
 * @details     Computes Euler/Navier-Stokes flux vector.  
 *              Jacobian matrix is replaced by first-order flux function,
 *              typically Rusanov flux is implemented.  
 *              Therefore, only convective flux is used.
 */
__device__ void ns_flux_container(UCFDReal *u, UCFDReal *nf, UCFDReal *f);

/**
 * @details     Computes convective flux for RANS one- or two-equations.  
 *              Similar to the Navier-Stokes equations,
 *              RANS equations can be reformulated into the finite-volume framework.  
 *              It contains conservative variables, convective/viscous flux, and source term.
 *              Convective flux in RANS equations is computed
 *              simply by multiplying conservative variables and contravariant velocity.
 */
__device__ void rans_flux_container(UCFDReal *u, UCFDReal *nf, UCFDReal *f);


ucfd_status_t rans_source_jacobian(UCFDReal *uf, UCFDReal tmat[NTURBVARS][NTURBVARS], UCFDReal *dsrc);

#endif