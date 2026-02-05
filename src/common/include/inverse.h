/**
 * @file        inverse.h
 * @brief       Header file for LU Decomposition/Substitution
 */
#ifndef INVERSE_H
#define INVERSE_H
#include "config.h"


/**
 * @brief       LU Decomposition function.
 * @param       n   Dimension of the A matrix
 * @param       A   Target matrix to decompose
 */
void ludcmp(int n, UCFD_FLOAT *A);

/**
 * @brief       Forward/Backward Substitution function.
 * @param       n   Dimension of the matrix
 * @param       LU  LU decomposed matrix
 * @param       b   Right-hand-side vector
 */
void lusub(int n, UCFD_FLOAT *LU, UCFD_FLOAT *b);

void lusubmat(int n, UCFD_FLOAT *LU, UCFD_FLOAT *B);

void lusubmattrans(int n, double *LU, double *B);


#endif //INVERSE_H