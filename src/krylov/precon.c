/** ======================================================================================================================
 * @file        precon.c
 * @brief       Preconditioners for Krylov subspace methods
 * @details     Preconditioner plays a crucial role in Krylov subspace methods, improving convergence characteristics.
 *              There are two preconditioner types, `BILU` and `LU-SGS`.
 *
 *              (1) BILU (Block fill-in Incomplete LU)
 *                  fill-in process is allowed in each block.
 *
 *              (2) LU-SGS (Lower-Upper Symmetric Gauss-Seidel)
 *                  Preconditioning version of LU-SGS method, modified for BSR format.
 * @author
 *              - Namhyoung Kim (knhkse@inha.edu), Department of Aerospace Engineering, Inha University
 *              - Jin Seok Park (jinseok.park@inha.ac.kr), Department of Aerospace Engineering, Inha University
 *
 * @date        April 2025
 * @version     1.0
 * @par         Copyright
 *              Copyright (c) 2024, Namhyoung Kim and Jin Seok Park, Inha University, All rights reserved.
 * @par         License
 *              This project is release under the terms of the MIT License (see LICENSE file).
 *
 * =======================================================================================================================
 */
#include "precon.h"
#include "inverse.h"
#include <stdlib.h>
#include <omp.h>

#define blkdim BLOCK*BLOCK

/**
 * @details     This function refactors non-zero values of BSR matrix applying block fill-in process.
 */
ucfd_status_t bilu_prepare(int bn, int *iw,
                           int *row_ptr, int *col_ind, int *diag_ind, double *nnz_data)
{
    int idx, kdx, ck, row, col, ele;
    int st, ed, jed, kk, kst, ked, jj, iwj;
    double v, Aik[BLOCK][BLOCK];

    for (idx = 0; idx < bn; idx++)
    {
        st = row_ptr[idx];
        ed = diag_ind[idx];
        jed = row_ptr[idx+1];

        // kdx : A[i,k]
        // kk : A[k,k]
        for (kdx = st; kdx < ed; kdx++)
        {
            ck = col_ind[kdx];
            kk = diag_ind[ck];
            kst = row_ptr[ck];
            ked = row_ptr[ck+1];
            
            // A[i,k] := A[i,k] @ inv(A[k,k])
            lusubmattrans(&nnz_data[kk*blkdim], &nnz_data[kdx*blkdim]);
            // memcpy(Aik, &nnz_data[kdx*blkdim], sizeof(double)*BLOCK);
            for (row=0; row<BLOCK; row++) {
                for (col=0; col<BLOCK; col++)
                    Aik[row][col] = nnz_data[kdx*blkdim+row*BLOCK+col];
            }

            // Prepare iw
            for (jj=kst; jj<ked; jj++) iw[col_ind[jj]] = jj;

            // j iteration
            for (jj=kdx+1; jj<jed; jj++) {
                iwj = iw[col_ind[jj]];

                if (iwj != 0) {
                    // nnz_data[jj] -= Aik * nnz_data[iwj]
                    for (row=0; row<BLOCK; row++) {
                        for (col=0; col<BLOCK; col++) {
                            v = 0.0;
                            for (ele=0; ele<BLOCK; ele++)
                                v += Aik[row][ele] * nnz_data[iwj*blkdim+ele*BLOCK+col];
                            nnz_data[jj*blkdim+row*BLOCK+col] -= v;
                        }
                    }
                }
            }

            // Clean iw
            for (jj=kst; jj<ked; jj++) iw[col_ind[jj]] = 0;
        }
        // Inverse current row diagonal matrix
        ludcmp(&nnz_data[ed*blkdim]);
    }

    return UCFD_STATUS_SUCCESS;
}


/**
 * @details     This function applies preconditioner matrix into arbitrary input vector `b`.
 *              In other words, solve `Px = b`.
 */
void bilu_psolve(int bn, int *row_ptr,
                 int *col_ind, int *diag_ind, double *nnz_data, double *b)
{
    int idx, jdx, kdx, row, col;
    int dd, st, ed, cind;
    double v, arr[BLOCK];

    // Forward substitution
    for (idx = 0; idx < bn; idx++)
    {
        dd = diag_ind[idx];
        st = row_ptr[idx];

        // Initialize arr
        for (kdx = 0; kdx < BLOCK; kdx++)
            arr[kdx] = b[kdx + idx * BLOCK];

        for (jdx = st; jdx < dd; jdx++)
        {
            cind = col_ind[jdx];

            for (row = 0; row < BLOCK; row++)
            {
                v = 0.0;
                for (col = 0; col < BLOCK; col++)
                    v += nnz_data[col + row * BLOCK + jdx * blkdim] * b[col + cind * BLOCK];
                arr[row] -= v;
            }
        }

        for (kdx = 0; kdx < BLOCK; kdx++)
            b[kdx + idx * BLOCK] = arr[kdx];
    }

    // Backward substitution
    for (idx = bn - 1; idx > -1; idx--)
    {
        dd = diag_ind[idx];
        ed = row_ptr[idx + 1];

        // Initialize
        for (kdx = 0; kdx < BLOCK; kdx++)
            arr[kdx] = b[kdx + idx * BLOCK];

        for (jdx = dd + 1; jdx < ed; jdx++)
        {
            cind = col_ind[jdx];

            for (row = 0; row < BLOCK; row++)
            {
                v = 0.0;
                for (col = 0; col < BLOCK; col++)
                    v += nnz_data[col + row * BLOCK + jdx * blkdim] * b[col + cind * BLOCK];
                arr[row] -= v;
            }
        }

        // LU substitution for vector
        lusub(&nnz_data[dd*blkdim], arr);
        for (row=0; row<BLOCK; row++) b[idx*BLOCK+row] = arr[row];
    }
}


/**
 * @details     LU decomposition is applied in every diagonal matrix.
 */
ucfd_status_t lusgs_prepare(int bn, int *diag_ind, double *nnz_data)
{
    int idx, didx;

    // Get diagonal block and store reverse
    // Parallel computation available (Only diagonal matrices are used)
    #pragma omp parallel for private(didx)
    for (idx = 0; idx < bn; idx++)
    {
        didx = diag_ind[idx];
        ludcmp(&nnz_data[didx * blkdim]);
    }

    return UCFD_STATUS_SUCCESS;
}

/**
 * @details     This function applies preconditioner matrix into arbitrary input vector `b`.
 *              In other words, solve `Px = b`.
 */
void lusgs_psolve(int bn, int *row_ptr,
                  int *col_ind, int *diag_ind, double *nnz_data, double *b)
{
    int idx, jdx, kdx, row, col;
    int dd, st, ed, cind;
    double v, arr[BLOCK];

    // Forward sweep : (D+L)x' = b -> x' = inv(D) * (b-Lx')
    for (idx = 0; idx < bn; idx++)
    {
        dd = diag_ind[idx];
        st = row_ptr[idx];

        // arr := b
        for (kdx = 0; kdx < BLOCK; kdx++)
            arr[kdx] = b[kdx + idx * BLOCK];

        // arr := b - Lx'
        for (jdx = st; jdx < dd; jdx++)
        {
            cind = col_ind[jdx];
            for (row = 0; row < BLOCK; row++)
            {
                v = 0.0;
                for (col = 0; col < BLOCK; col++)
                    v += nnz_data[col + row * BLOCK + jdx * blkdim] * b[col + cind * BLOCK];
                arr[row] -= v;
            }
        }

        // x' := inv(D) * (b-Lx') = inv(D) * arr
        lusub(&nnz_data[dd * blkdim], arr);
        for (kdx = 0; kdx < BLOCK; kdx++)
            b[kdx + idx * BLOCK] = arr[kdx];
    }

    // Backward sweep : (D+U)x = Dx' -> x = x' - inv(D) * Ux
    for (idx = bn - 1; idx > -1; idx--)
    {
        dd = diag_ind[idx];
        ed = row_ptr[idx + 1];

        // Initialize
        for (kdx = 0; kdx < BLOCK; kdx++)
            arr[kdx] = 0.0;

        // arr := Ux
        for (jdx = dd + 1; jdx < ed; jdx++)
        {
            cind = col_ind[jdx];
            for (row = 0; row < BLOCK; row++)
            {
                v = 0.0;
                for (col = 0; col < BLOCK; col++)
                    v += nnz_data[col + row * BLOCK + jdx * blkdim] * b[col + cind * BLOCK];
                arr[row] += v;
            }
        }

        // arr := inv(D) Ux
        lusub(&nnz_data[dd * blkdim], arr);

        // b := b - inv(D) Ux
        for (kdx = 0; kdx < BLOCK; kdx++)
            b[kdx + idx * BLOCK] -= arr[kdx];
    }
}

void none_psolve(int bn, int *row_ptr,
                 int *col_ind, int *diag_ind, double *nnz_data, double *b) {};