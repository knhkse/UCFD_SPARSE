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
void ludcmp(UCFDInt block, UCFDReal *A)
{
    UCFDInt row, col, kdx, nrow;
    UCFDReal val;

    if (block == 1) {               // 1-equation RANS model
        A[0] = 1.0/A[0];
    }

    else {
        for (row=1; row<block; row++) {     // Skip first row
            nrow = block*row;
            A[nrow] /= A[0];
            for (col=1; col<block; col++) {
                // Lower triangular matrix
                if (row > col) {
                    val = 0.0;
                    for (kdx=0; kdx<col; kdx++)
                        val += A[nrow+kdx] * A[col + block*kdx];
                    A[nrow+col] = (A[nrow+col] - val)/A[(block+1)*col];
                }

                // Upper triangular matrix
                else {
                    val = 0.0;
                    for (kdx=0; kdx<row; kdx++)
                        val += A[nrow+kdx]*A[block*kdx+col];
                    A[nrow+col] -= val;
                }
            }
        }
    }
}


/**
 * @details     This function performs Forward/Backward substitution of LU decomposed matrix.
 */
void lusub(UCFDInt block, UCFDReal *LU, UCFDReal *b)
{
    UCFDInt row, col, nrow;
    UCFDReal val;

    if (block == 1) {                       // 1-equation RANS model
        b[0] *= LU[0];
    }

    else {
        // Forward substitution
        for (row=1; row<block; row++) {
            nrow = row*block;
            val = 0.0;
            for (col=0; col<row; col++)
                val += LU[nrow+col]*b[col];
            b[row] -= val;
        }

        // Backward substitution
        b[block-1] /= LU[block*block-1];
        for (row=block-2; row>-1; row--) {
            nrow = row*block;
            val = 0.0;
            for (col=row+1; col<block; col++)
                val += LU[nrow+col]*b[col];
            b[row] = (b[row] - val)/LU[nrow+row];
        }
    }
}

void lusubmattrans(UCFDInt block, UCFDReal *LU, UCFDReal *B)
{
    UCFDInt row, col, scol;
    UCFDReal val;

    if (block == 1) {                       // 1-equation RANS model
        B[0] *= LU[0];
    }

    else {
        // Forward substitution
        for (scol=0; scol<block; scol++) B[scol*block] /= LU[0];
        for (row=1; row<block; row++) {
            for (scol=0; scol<block; scol++) {
                val = 0.0;
                for (col=0; col<row; col++)
                    val += B[scol*block+col] * LU[col*block+row];
                B[scol*block+row] = (B[scol*block+row] - val)/LU[row*block+row];
            }
        }

        // Backward substitution
        for (row=block-2; row>-1; row--) {
            for (scol=0; scol<block; scol++) {
                val = 0.0;
                for (col=row+1; col<block; col++)
                    val += B[scol*block+col] * LU[col*block+row];
                B[scol*block+row] -= val;
            }
        }
    }
}
