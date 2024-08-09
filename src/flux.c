/** ======================================================================================================================
 * @file        flux.c
 * @brief       Computes numerical flux from conservative variables.  
 *              Non-physical value correction included.
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

#include "flux.h"

/** Default value */
#ifndef CONST_H
    #define gamma   1.4     /* Specific heat ratio */
    #define pmin    1e-15   /* Minimum pressure value */
#endif


/**
 * @details     Computes Euler/Navier-Stokes flux vector.  
 *              Jacobian matrix is replaced by first-order flux function,
 *              typically Rusanov flux is implemented.  
 *              Therefore, only convective flux is used.
 */
void ns_flux_container(int nfvars, int ndims, double *u, double *nf, double *f)
{
    double rho = u[0];          /** density */
    double et = u[nfvars-1];    /** Total Energy */
    double temp = 0.0;          /** \f$\rho^2 \times (u^2 + v^2)\f$ */
    double contrav = 0.0;       /** Contravariant velocity */
    int i;

    for (i=0; i<ndims; i++) {
        contrav += u[i+1]*nf[i];
        temp += u[i+1]*u[i+1];
    }
    contrav /= rho;

    // Apply lower bound of pressure value
    double p = (gamma - 1.0)*(et - 0.5*temp/rho);
    if (p < pmin) {
        p = pmin;
        et = p/(gamma-1.0) + 0.5*temp/rho;
        u[nfvars-1] = et;
    }
    
    // Total enthalpy
    double ht = et + p;

    // Computes flux array
    f[0] = rho*contrav;
    for (int i=0; i<ndims; i++) {
        f[i+1] = u[i+1] * contrav + nf[i]*p;
    }
    f[nfvars-1] = ht*contrav;
}

/**
 * @details     Computes convective flux for RANS one- or two-equations.  
 *              Similar to the Navier-Stokes equations,
 *              RANS equations can be reformulated into the finite-volume framework.  
 *              It contains conservative variables, convective/viscous flux, and source term.
 *              Convective flux in RANS equations is computed
 *              simply by multiplying conservative variables and contravariant velocity.
 */
void rans_flux_container(int nfvars, int ndims, int nturbvars, double *u, double *nf, double *f)
{
    double rho = u[0];
    double contrav = 0.0;

    for (int i=0; i<ndims; i++) {
        contrav += u[i+1] * nf[i];
    }
    contrav /= rho;

    for (int i=0; i<nturbvars; i++) {
        f[i] = u[nfvars+i]*contrav;
    }
}