/**
 * @file        inverse.h
 * @brief       Header file for LU decomposition/substitution
 */
#pragma once

#include "config.h"


/**
 * @brief       LU Decomposition function
 * @param       A   Target matrix to decompose
 */
void ludcmp(UCFDInt block, UCFDReal *A);

/**
 * @brief       Forward/Backward substitution function
 * @param       LU  LU decomposed matrix
 * @param       b   Right-hand-side vector
 */
void lusub(UCFDInt block, UCFDReal *LU, UCFDReal *b);

/**
 * @brief       Forward/Backward substitution for transposed matrix
 * @param       LU  LU decomposed matrix
 * @param       B   Right-hand side matrix
 */
void lusubmattrans(UCFDInt block, UCFDReal *LU, UCFDReal *B);
