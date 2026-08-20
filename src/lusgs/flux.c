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
void ns_flux_container(UCFDInt nfvars, UCFDInt nturbvars, UCFDInt ndims, UCFDReal *u, UCFDReal *nf, UCFDReal *f)
{   
    /**
     * Variable description :  
     * `rho` : Density  
     * `et` : Total Energy  
     * `temp` : \f$\rho^2 \times (u^2 + v^2)\f$  
     * `contrav` : Contravariant velocity
     */
    UCFDReal rho = u[0];
    UCFDReal et = u[nfvars-1];
    UCFDReal temp = 0.0;
    UCFDReal contrav = 0.0;
    int i;

    for (i=0; i<ndims; ++i) {
        contrav += u[i+1]*nf[i];
        temp += u[i+1]*u[i+1];
    }
    contrav /= rho;

    // Apply lower bound of pressure value
    UCFDReal p = (GAMMA - 1.0)*(et - 0.5*temp/rho);
    if (p < PMIN) {
        p = PMIN;
        et = p/(GAMMA-1.0) + 0.5*temp/rho;
        u[nfvars-1] = et;
    }
    
    // Total enthalpy
    UCFDReal ht = et + p;

    // Computes flux array
    f[0] = rho*contrav;
    for (UCFDInt i=0; i<ndims; ++i) {
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
void rans_flux_container(UCFDInt nfvars, UCFDInt nturbvars, UCFDInt ndims, UCFDReal *u, UCFDReal *nf, UCFDReal *f)
{
    UCFDReal rho = u[0];
    UCFDReal contrav = 0.0;

    for (UCFDInt i=0; i<ndims; ++i) {
        contrav += u[i+1] * nf[i];
    }
    contrav /= rho;

    for (UCFDInt i=0; i<nturbvars; ++i) {
        f[i] = u[nfvars+i]*contrav;
    }
}

