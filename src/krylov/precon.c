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

#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "mkl.h"
#include "precon.h"
#include "inverse.h"

/**
 * @details     This function refactors non-zero values of BSR matrix applying block fill-in process.
 */
ucfd_status_t bilu_prepare(const int bn, const int blk, int *iw,
                           const int *row_ptr, const int *col_ind, const int *diag_ind, double *nnz_data)
{
    int idx, jdx, kdx, ck, mdx, row, col, ele;
    int st, ed, jed, kk, kst, ked, jj, iwj;
    double v, Aik[blk * blk];
    int blk2 = blk * blk;

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
            lumatsubtrans(blk, &nnz_data[kk*blk2], &nnz_data[kdx*blk2]);
            // memcpy(Aik, &nnz_data[kdx*blk2], sizeof(double)*blk);
            for (row=0; row<blk; row++) {
                for (col=0; col<blk; col++)
                    Aik[row*blk+col] = nnz_data[kdx*blk2+row*blk+col];
            }

            // Prepare iw
            for (jj=kst; jj<ked; jj++) iw[col_ind[jj]] = jj;

            // j iteration
            for (jj=kdx+1; jj<jed; jj++) {
                iwj = iw[col_ind[jj]];

                if (iwj != 0) {
                    // nnz_data[jj] -= Aik * nnz_data[iwj]
                    for (row=0; row<blk; row++) {
                        for (col=0; col<blk; col++) {
                            v = 0.0;
                            for (ele=0; ele<blk; ele++)
                                v += Aik[row*blk+ele] * nnz_data[iwj*blk2+ele*blk+col];
                            nnz_data[jj*blk2+row*blk+col] -= v;
                        }
                    }
                }
            }

            // Clean iw
            for (jj=kst; jj<ked; jj++) iw[col_ind[jj]] = 0;
        }
        // Inverse current row diagonal matrix
        ludcmp(blk, &nnz_data[ed*blk2]);
    }

    return UCFD_STATUS_SUCCESS;
}


/**
 * @details     This function applies preconditioner matrix into arbitrary input vector `b`.
 *              In other words, solve `Px = b`.
 */
ucfd_status_t bilu_psolve(const int bn, const int blk, const int *row_ptr,
                          const int *col_ind, const int *diag_ind, double *nnz_data, double *b)
{
    int idx, jdx, kdx, row, col;
    int dd, st, ed, cind;
    double v, arr[blk];
    const int blk2 = blk*blk;

    // Forward substitution
    for (idx = 0; idx < bn; idx++)
    {
        dd = diag_ind[idx];
        st = row_ptr[idx];

        // Initialize arr
        for (kdx = 0; kdx < blk; kdx++)
            arr[kdx] = b[kdx + idx * blk];

        for (jdx = st; jdx < dd; jdx++)
        {
            cind = col_ind[jdx];

            for (row = 0; row < blk; row++)
            {
                v = 0.0;
                for (col = 0; col < blk; col++)
                    v += nnz_data[col + row * blk + jdx * blk2] * b[col + cind * blk];
                arr[row] -= v;
            }
        }

        for (kdx = 0; kdx < blk; kdx++)
            b[kdx + idx * blk] = arr[kdx];
    }

    // Backward substitution
    for (idx = bn - 1; idx > -1; idx--)
    {
        dd = diag_ind[idx];
        ed = row_ptr[idx + 1];

        // Initialize
        for (kdx = 0; kdx < blk; kdx++)
            arr[kdx] = b[kdx + idx * blk];

        for (jdx = dd + 1; jdx < ed; jdx++)
        {
            cind = col_ind[jdx];

            for (row = 0; row < blk; row++)
            {
                v = 0.0;
                for (col = 0; col < blk; col++)
                    v += nnz_data[col + row * blk + jdx * blk2] * b[col + cind * blk];
                arr[row] -= v;
            }
        }

        // LU substitution for vector
        luvecsub(blk, &nnz_data[dd*blk2], arr);
        for (row=0; row<blk; row++) b[idx*blk+row] = arr[row];
    }

    return UCFD_STATUS_SUCCESS;
}


/**
 * @details     LU decomposition is applied in every diagonal matrix.
 */
ucfd_status_t lusgs_prepare(const int bn, const int blk, const int *diag_ind, double *nnz_data)
{
    const int blk2 = blk * blk;
    int idx, didx;

    // Get diagonal block and store reverse
    // Parallel computation available (Only diagonal matrices are used)
    #pragma omp parall for private(didx)
    for (idx = 0; idx < bn; idx++)
    {
        didx = diag_ind[idx];
        ludcmp(blk, &nnz_data[didx * blk2]);
    }

    return UCFD_STATUS_SUCCESS;
}

/**
 * @details     This function applies preconditioner matrix into arbitrary input vector `b`.
 *              In other words, solve `Px = b`.
 */
ucfd_status_t lusgs_psolve(const int bn, const int blk, const int *row_ptr,
                           const int *col_ind, const int *diag_ind, double *nnz_data, double *b)
{
    const int blk2 = blk * blk;
    int idx, jdx, kdx, row, col;
    int dd, st, ed, cind;
    double v, arr[blk];

    // Forward sweep : (D+L)x' = b -> x' = inv(D) * (b-Lx')
    for (idx = 0; idx < bn; idx++)
    {
        dd = diag_ind[idx];
        st = row_ptr[idx];

        // arr := b
        for (kdx = 0; kdx < blk; kdx++)
            arr[kdx] = b[kdx + idx * blk];

        // arr := b - Lx'
        for (jdx = st; jdx < dd; jdx++)
        {
            cind = col_ind[jdx];
            for (row = 0; row < blk; row++)
            {
                v = 0.0;
                for (col = 0; col < blk; col++)
                    v += nnz_data[col + row * blk + jdx * blk2] * b[col + cind * blk];
                arr[row] -= v;
            }
        }

        // x' := inv(D) * (b-Lx') = inv(D) * arr
        lusubst(blk, &nnz_data[dd * blk2], arr);
        for (kdx = 0; kdx < blk; kdx++)
            b[kdx + idx * blk] = arr[kdx];
    }

    // Backward sweep : (D+U)x = Dx' -> x = x' - inv(D) * Ux
    for (idx = bn - 1; idx > -1; idx--)
    {
        dd = diag_ind[idx];
        ed = row_ptr[idx + 1];

        // Initialize
        for (kdx = 0; kdx < blk; kdx++)
            arr[kdx] = 0.0;

        // arr := Ux
        for (jdx = dd + 1; jdx < ed; jdx++)
        {
            cind = col_ind[jdx];
            for (row = 0; row < blk; row++)
            {
                v = 0.0;
                for (col = 0; col < blk; col++)
                    v += nnz_data[col + row * blk + jdx * blk2] * b[col + cind * blk];
                arr[row] += v;
            }
        }

        // arr := inv(D) Ux
        lusubst(blk, &nnz_data[dd * blk2], arr);

        // b := b - inv(D) Ux
        for (kdx = 0; kdx < blk; kdx++)
            b[kdx + idx * blk] -= arr[kdx];
    }

    return UCFD_STATUS_SUCCESS;
}
