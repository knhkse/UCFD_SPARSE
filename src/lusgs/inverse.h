#ifndef INVERSE_H
#define INVERSE_H

/**
 * @file        inverse.h
 * @brief       Header file for LU Decomposition/Substitution
 */

/**
 * @brief       LU Decomposition function.
 * @param       n   Dimension of the A matrix
 * @param       A   Target matrix to decompose
 */
void ludcmp(int n, double *A);

/**
 * @brief       Forward/Backward Substitution function.
 * @param       n   Dimension of the matrix
 * @param       LU  LU decomposed matrix
 * @param       b   Right-hand-side vector
 */
void lusubst(int n, double *LU, double *b);


#endif //INVERSE_H