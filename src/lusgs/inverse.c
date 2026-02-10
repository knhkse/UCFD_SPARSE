/** ======================================================================================================================
 * @file        inverse.c
 * @brief       Compute matrix inverse using `LU Decomposition`, and forward/backward substitution.
 * @details     This file contains LU Decomposition function and substitution function.  
 *                
 *              (1) LU Decomposition : ludcmp  
 *              Decompose the input matrix A by Lower/Upper triangular matrix.  
 *              A = LU  
 *                
 *              (2) Substitution : lusubst  
 *              Solve Ax = b by using forward/backward substitution.  
 *              Input vector `b` is overwritten with the solution vector `x`.  
 * 
 * @note        Input matrix must be the n-by-n square matrix, and `Row-major` format.  
 *              Each function treats target matrix as a one-dimensional array.
 * 
 * @author
 *              - Namhyoung Kim (knhkse@inha.edu), Department of Aerospace Engineering, Inha University
 *              - Jin Seok Park (jinseok.park@inha.ac.kr), Department of Aerospace Engineering, Inha University
 * 
 * @date        Nov 2024
 * @version     1.0
 * @par         Copyright
 *              Copyright (c) 2024, Namhyoung Kim and Jin Seok Park, Inha University, All rights reserved.
 * @par         License
 *              This project is release under the terms of the MIT License (see LICENSE file).
 * 
 * =======================================================================================================================
 */
#include "inverse.h"

// TODO : Bespoke code generation

/**
 * @details     Decompose matrix A into lower and upper triangular matrix
 */
void ludcmp(int n, UCFD_FLOAT A[n][n])
{
    int row, col, kdx;
    UCFD_FLOAT val;

    if (n == 1) {               // 1-equation RANS model
        A[0][0] = 1.0/A[0][0];
    }

    else {
        for (row=1; row<n; row++) {     // Skip first row
            A[row][0] /= A[0][0];
            for (col=1; col<n; col++) {
                // Lower triangular matrix
                if (row > col) {
                    val = 0.0;
                    for (kdx=0; kdx<col; kdx++)
                        val += A[row][kdx] * A[kdx][col];
                    A[row][col] = (A[row][col] - val)/A[col][col];
                }

                // Upper triangular matrix
                else {
                    val = 0.0;
                    for (kdx=0; kdx<row; kdx++)
                        val += A[row][kdx]*A[kdx][col];
                    A[row][col] -= val;
                }
            }
        }
    }
}


/**
 * @details     This function performs Forward/Backward substitution of LU decomposed matrix.
 */
void lusub(int n, UCFD_FLOAT LU[n][n], UCFD_FLOAT *b)
{
    int row, col;
    UCFD_FLOAT val;

    if (n == 1) {                       // 1-equation RANS model
        b[0] *= LU[0][0];
    }

    else {
        // Forward substitution
        for (row=1; row<n; row++) {
            val = 0.0;
            for (col=0; col<row; col++)
                val += LU[row][col]*b[col];
            b[row] -= val;
        }

        // Backward substitution
        b[n-1] /= LU[n-1][n-1];
        for (row=n-2; row>-1; row--) {
            val = 0.0;
            for (col=row+1; col<n; col++)
                val += LU[row][col]*b[col];
            b[row] = (b[row] - val)/LU[row][row];
        }
    }
}