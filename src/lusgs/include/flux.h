/**
 * @file        flux.h
 * @brief       Header file for numerical flux funtions
 * @details     Declaration of convective flux for Navier-Stokes and RANS equations.
 */
#ifndef FLUX_H
#define FLUX_H
#include "ucfd_types.h"
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
 * @brief       Computes flux for Navier-Stokes equations.
 * @param       u           Conservative vector
 * @param       nf          Surface vector
 * @param       f           Flux vector
 */
void ns_flux_container(UCFD_FLOAT *u, UCFD_FLOAT *nf, UCFD_FLOAT *f);

/**
 * @brief       Computes flux for RANS equations.
 * @param       u           Conservative vector
 * @param       nf          Surface vector
 * @param       f           Flux vector
 */
void rans_flux_container(UCFD_FLOAT *u, UCFD_FLOAT *nf, UCFD_FLOAT *f);


/**
 * @brief       Computes source term Jacobian matrix for RANS equations.
 * @param       betast      beta* value for kw-SST RANS model
 * @param       uf          Conservative vector
 * @param       tmat        Turbulence Jacobian matrix
 * @param       dsrc        Source term derivatives vector
 */
ucfd_status_t rans_source_jacobian(UCFD_FLOAT *uf, UCFD_FLOAT *tmat[], UCFD_FLOAT *dsrc);

void print_configure(); // Remove later

#endif //FLUX_H