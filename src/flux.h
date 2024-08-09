#ifndef FLUX_H
#define FLUX_H

/**
 * @brief       Computes flux for Navier-Stokes equations.
 * @param       nfvars      Number of flux variables
 * @param       ndims       Dimensions
 * @param       u           Conservative vector
 * @param       nf          Surface vector
 * @param       f           Flux vector
 */
void ns_flux_container(int nfvars, int ndims, double *u, double *nf, double *f);


/**
 * @brief       Computes flux for RANS equations.
 * @param       nfvars      Number of flux variables
 * @param       ndims       Dimensions
 * @param       nturbvars   Number of turbulence variables
 * @param       u           Conservative vector
 * @param       nf          Surface vector
 * @param       f           Flux vector
 */
void rans_flux_container(int nfvars, int ndims, int nturbvars, double *u, double *nf, double *f);


#endif //FLUX_H