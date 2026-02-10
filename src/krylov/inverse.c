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
void ludcmp(UCFD_FLOAT *A)
{
    int row, col, kdx, nrow;
    UCFD_FLOAT val;

    if (BLOCK == 1) {               // 1-equation RANS model
        A[0] = 1.0/A[0];
    }

    else {
        for (row=1; row<BLOCK; row++) {     // Skip first row
            nrow = BLOCK*row;
            A[nrow] /= A[0];
            for (col=1; col<BLOCK; col++) {
                // Lower triangular matrix
                if (row > col) {
                    val = 0.0;
                    for (kdx=0; kdx<col; kdx++)
                        val += A[nrow+kdx] * A[col + BLOCK*kdx];
                    A[nrow+col] = (A[nrow+col] - val)/A[(BLOCK+1)*col];
                }

                // Upper triangular matrix
                else {
                    val = 0.0;
                    for (kdx=0; kdx<row; kdx++)
                        val += A[nrow+kdx]*A[BLOCK*kdx+col];
                    A[nrow+col] -= val;
                }
            }
        }
    }
}


/**
 * @details     This function performs Forward/Backward substitution of LU decomposed matrix.
 */
void lusub(UCFD_FLOAT *LU, UCFD_FLOAT *b)
{
    int row, col, nrow;
    UCFD_FLOAT val;

    if (BLOCK == 1) {                       // 1-equation RANS model
        b[0] *= LU[0];
    }

    else {
        // Forward substitution
        for (row=1; row<BLOCK; row++) {
            nrow = row*BLOCK;
            val = 0.0;
            for (col=0; col<row; col++)
                val += LU[nrow+col]*b[col];
            b[row] -= val;
        }

        // Backward substitution
        b[BLOCK-1] /= LU[BLOCK*BLOCK-1];
        for (row=BLOCK-2; row>-1; row--) {
            val = 0.0;
            for (col=row+1; col<BLOCK; col++)
                val += LU[nrow+col]*b[col];
            b[row] = (b[row] - val)/LU[nrow+row];
        }
    }
}

void lusubmattrans(UCFD_FLOAT *LU, UCFD_FLOAT *B)
{
    int row, col, scol;
    UCFD_FLOAT val;

    if (BLOCK == 1) {                       // 1-equation RANS model
        B[0] *= LU[0];
    }

    else {
        // Forward substitution
        for (scol=0; scol<BLOCK; scol++) B[scol*BLOCK] /= LU[0];
        for (row=1; row<BLOCK; row++) {
            for (scol=0; scol<BLOCK; scol++) {
                val = 0.0;
                for (col=0; col<row; col++)
                    val += B[scol*BLOCK+col] * LU[col*BLOCK+row];
                B[scol*BLOCK+row] = (B[scol*BLOCK+row] - val)/LU[row*BLOCK+row];
            }
        }

        // Backward substitution
        for (row=BLOCK-2; row>-1; row--) {
            for (scol=0; scol<BLOCK; scol++) {
                val = 0.0;
                for (col=row+1; col<BLOCK; col++)
                    val += B[scol*BLOCK+col] * LU[col*BLOCK+row];
                B[scol*BLOCK+row] -= val;
            }
        }
    }
}