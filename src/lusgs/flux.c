/** ======================================================================================================================
 * @file        flux.c
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
#include "flux.h"

#define max(a,b) (((a) > (b)) ? (a) : (b))

/**
 * @details     Computes Euler/Navier-Stokes flux vector.  
 *              Jacobian matrix is replaced by first-order flux function,
 *              typically Rusanov flux is implemented.  
 *              Therefore, only convective flux is used.
 */
void ns_flux_container(UCFD_FLOAT *u, UCFD_FLOAT *nf, UCFD_FLOAT *f)
{   
    /**
     * Variable description :  
     * `rho` : Density  
     * `et` : Total Energy  
     * `temp` : \f$\rho^2 \times (u^2 + v^2)\f$  
     * `contrav` : Contravariant velocity
     */
    UCFD_FLOAT rho = u[0];
    UCFD_FLOAT et = u[NFVARS-1];
    UCFD_FLOAT temp = 0.0;
    UCFD_FLOAT contrav = 0.0;
    int i;

    for (i=0; i<NDIMS; i++) {
        contrav += u[i+1]*nf[i];
        temp += u[i+1]*u[i+1];
    }
    contrav /= rho;

    // Apply lower bound of pressure value
    UCFD_FLOAT p = (GAMMA - 1.0)*(et - 0.5*temp/rho);
    if (p < PMIN) {
        p = PMIN;
        et = p/(GAMMA-1.0) + 0.5*temp/rho;
        u[NFVARS-1] = et;
    }
    
    // Total enthalpy
    UCFD_FLOAT ht = et + p;

    // Computes flux array
    f[0] = rho*contrav;
    for (UCFD_INT i=0; i<NDIMS; i++) {
        f[i+1] = u[i+1] * contrav + nf[i]*p;
    }
    f[NFVARS-1] = ht*contrav;
}

/**
 * @details     Computes convective flux for RANS one- or two-equations.  
 *              Similar to the Navier-Stokes equations,
 *              RANS equations can be reformulated into the finite-volume framework.  
 *              It contains conservative variables, convective/viscous flux, and source term.
 *              Convective flux in RANS equations is computed
 *              simply by multiplying conservative variables and contravariant velocity.
 */
void rans_flux_container(UCFD_FLOAT *u, UCFD_FLOAT *nf, UCFD_FLOAT *f)
{
    UCFD_FLOAT rho = u[0];
    UCFD_FLOAT contrav = 0.0;

    for (UCFD_INT i=0; i<NDIMS; i++) {
        contrav += u[i+1] * nf[i];
    }
    contrav /= rho;

    for (UCFD_INT i=0; i<NTURBVARS; i++) {
        f[i] = u[NFVARS+i]*contrav;
    }
}


ucfd_status_t rans_source_jacobian(UCFD_FLOAT *uf, UCFD_FLOAT tmat[NTURBVARS][NTURBVARS], UCFD_FLOAT *dsrc)
{
    /* 1-equation RANS model (Spalart-Allmaras) */
    if (NTURBVARS == 1) tmat[0][0] += dsrc[NVARS-1];

    /* 2-equations RANS model (kw-SST) */
    else if (NTURBVARS == 2) {
        UCFD_FLOAT k = uf[NVARS-2] / uf[0];
        tmat[0][0] += dsrc[NVARS-2];
        tmat[0][1] += max(BETAST*k, 0.0);
        tmat[1][1] += dsrc[NVARS-1];
    }
    else return UCFD_STATUS_NOT_SUPPORTED;

    return UCFD_STATUS_SUCCESS;
}