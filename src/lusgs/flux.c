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

// #ifndef CONST_H
//     /** Specific heat ratio */
//     #define gamma   1.4
//     /** Minimum pressure value */
//     #define pmin    1e-15
// #endif

#define max(a,b) (((a) > (b)) ? (a) : (b))

/**
 * @details     Computes Euler/Navier-Stokes flux vector.  
 *              Jacobian matrix is replaced by first-order flux function,
 *              typically Rusanov flux is implemented.  
 *              Therefore, only convective flux is used.
 */
void ns_flux_container(int nfvars, int ndims, double *u, double *nf, double *f)
{   
    /**
     * Variable description :  
     * `rho` : Density  
     * `et` : Total Energy  
     * `temp` : \f$\rho^2 \times (u^2 + v^2)\f$  
     * `contrav` : Contravariant velocity
     */
    double rho = u[0];
    double et = u[nfvars-1];
    double temp = 0.0;
    double contrav = 0.0;
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


int rans_source_jacobian(int nvars, int ntvars, double betast, \
                         double *uf, double *tmat, double *dsrc)
{
    /* 1-equation RANS model (Spalart-Allmaras) */
    if (ntvars == 1) {
        tmat[0] += dsrc[nvars-1];
    }

    /* 2-equations RANS model (kw-SST) */
    else if (ntvars == 2) {
        double k = uf[nvars-2] / uf[0];
        tmat[0] += dsrc[nvars-2];
        tmat[1] += max(betast*k, 0.0);
        tmat[3] += dsrc[nvars-1];
    }

    else {
        return -1;
    }

    return 0;
}