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
void ludcmp(int n, UCFD_FLOAT A[n][n]);

/**
 * @brief       Forward/Backward Substitution function.
 * @param       n   Dimension of the matrix
 * @param       LU  LU decomposed matrix
 * @param       b   Right-hand-side vector
 */
void lusub(int n, UCFD_FLOAT LU[n][n], UCFD_FLOAT *b);

#endif //INVERSE_H