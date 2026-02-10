/**
 * @file        inverse.h
 * @brief       Header file for LU decomposition/substitution
 */
#ifndef INVERSE_H
#define INVERSE_H
#include "config.h"


/**
 * @brief       LU Decomposition function
 * @param       A   Target matrix to decompose
 */
void ludcmp(UCFD_FLOAT *A);

/**
 * @brief       Forward/Backward substitution function
 * @param       LU  LU decomposed matrix
 * @param       b   Right-hand-side vector
 */
void lusub(UCFD_FLOAT *LU, UCFD_FLOAT *b);

/**
 * @brief       Forward/Backward substitution for transposed matrix
 * @param       LU  LU decomposed matrix
 * @param       B   Right-hand side matrix
 */
void lusubmattrans(UCFD_FLOAT *LU, UCFD_FLOAT *B);


#endif //INVERSE_H