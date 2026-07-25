/**
 * @file        inverse.h
 * @brief       Header file for LU decomposition/substitution
 */
#pragma once

#include "config.h"
#include "macros.h"

/**
 * @details     Decompose matrix A into lower and upper triangular matrix
 */
static void MAYBE_UNUSED ludcmp(UCFDInt block, UCFDReal *A)
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
static void MAYBE_UNUSED lusub(UCFDInt block, UCFDReal *LU, UCFDReal *b)
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

static void MAYBE_UNUSED lusubmattrans(UCFDInt block, UCFDReal *LU, UCFDReal *B)
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
