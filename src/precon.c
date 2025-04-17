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
#include "mkl.h"
#include "precon.h"
#include "inverse.h"

/**
 * @details     This function refactors non-zero values of BSR matrix applying block fill-in process.
 */
ucfd_status_t bilu_prepare(const int neles, const int nvars,
                           const int *row_ptr, const int *col_ind, const int *diag_ind, double *nnz_data)
{
    int idx, jdx, kdx, ldx, mdx, row, col, ele;
    int st, ed, ijdx, kjdx, ijed, kjed;
    double v, mat[nvars * nvars];
    int ijarr[nvars + 3], kjarr[nvars + 3];
    int nvars2 = nvars * nvars;
    lapack_int info, *ipiv;

    // Allocate memory
    ipiv = (lapack_int *)malloc(sizeof(lapack_int) * nvars);

    for (idx = 0; idx < neles; idx++) {
        st = row_ptr[idx];
        ed = diag_ind[idx];

        for (kdx = st; kdx < ed; kdx++) {
            ldx = diag_ind[col_ind[kdx]];

            for (row = 0; row < nvars; row++) {
                for (col = 0; col < nvars; col++) {
                    v = 0.0;
                    for (ele = 0; ele < nvars; ele++)
                        v += nnz_data[ele + row * nvars + kdx * nvars2] * nnz_data[col + ele * nvars + ldx * nvars2];
                    mat[col + row * nvars] = v;
                }
            }

            for (row = 0; row < nvars; row++) {
                for (col = 0; col < nvars; col++)
                    nnz_data[col + row * nvars + kdx * nvars2] = mat[col + row * nvars];
            }

            ijdx = kdx + 1;
            kjdx = ldx + 1;
            ijed = row_ptr[idx + 1];
            kjed = row_ptr[col_ind[kdx] + 1];
            mdx = 0;
            while (kjdx < kjed && ijdx < ijed) {
                if (col_ind[ijdx] == col_ind[kjdx]) {
                    ijarr[mdx] = ijdx;
                    kjarr[mdx] = kjdx;
                    mdx++;
                    kjdx++;
                    ijdx++;
                }
                else if (col_ind[ijdx] < col_ind[kjdx])
                    ijdx++;
                else
                    kjdx++;
            }

            // A[i,j] -= A[i,k] @ A[k,j]
            for (jdx = 0; jdx < mdx; jdx++) {
                ijdx = ijarr[jdx];
                kjdx = kjarr[jdx];
                for (row = 0; row < nvars; row++) {
                    for (col = 0; col < nvars; col++) {
                        v = 0.0;
                        for (ele = 0; ele < nvars; ele++)
                            v += mat[ele + row * nvars] * nnz_data[col + ele * nvars + kjdx * nvars2];
                        nnz_data[col + row * nvars + ijdx * nvars2] -= v;
                    }
                }
            }
        }
        // Inverse current row diagonal matrix
        info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, nvars, nvars, &nnz_data[ed * nvars * nvars], nvars, ipiv);
        if (info != 0)
            printf("Lapack dgetrf function error\n");
        info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, nvars, &nnz_data[ed * nvars * nvars], nvars, ipiv);
        if (info != 0)
            printf("Lapack dgetrf function error\n");
    }
    free(ipiv);

    return UCFD_STATUS_SUCCESS;
}


/**
 * @details     This function applies preconditioner matrix into arbitrary input vector `b`.
 *              In other words, solve `Px = b`.
 */
ucfd_status_t bilu_psolve(const int neles, const int nvars, const int *row_ptr,
                          const int *col_ind, const int *diag_ind, double *nnz_data, double *b)
{
    const int nvars2 = nvars * nvars;
    int idx, jdx, kdx, row, col;
    int dd, st, ed, cind;
    double v, arr[nvars];

    for (idx = 0; idx < neles; idx++)
    {
        dd = diag_ind[idx];
        st = row_ptr[idx];

        // Initialize arr
        for (kdx = 0; kdx < nvars; kdx++)
            arr[kdx] = b[kdx + idx * nvars];

        for (jdx = st; jdx < dd; jdx++)
        {
            cind = col_ind[jdx];

            for (row = 0; row < nvars; row++)
            {
                v = 0.0;
                for (col = 0; col < nvars; col++)
                    v += nnz_data[col + row * nvars + jdx * nvars2] * b[col + cind * nvars];
                arr[row] -= v;
            }
        }

        for (kdx = 0; kdx < nvars; kdx++)
            b[kdx + idx * nvars] = arr[kdx];
    }

    for (idx = neles - 1; idx > -1; idx--)
    {
        dd = diag_ind[idx];
        ed = row_ptr[idx + 1];

        // Initialize
        for (kdx = 0; kdx < nvars; kdx++)
            arr[kdx] = b[kdx + idx * nvars];

        for (jdx = dd + 1; jdx < ed; jdx++)
        {
            cind = col_ind[jdx];

            for (row = 0; row < nvars; row++)
            {
                v = 0.0;
                for (col = 0; col < nvars; col++)
                    v += nnz_data[col + row * nvars + jdx * nvars2] * b[col + cind * nvars];
                arr[row] -= v;
            }
        }

        // Multiply diagonal inverse
        for (row = 0; row < nvars; row++)
        {
            v = 0.0;
            for (col = 0; col < nvars; col++)
                v += nnz_data[col + row * nvars + dd * nvars2] * arr[col];
            b[row + idx * nvars] = v;
        }
    }

    return UCFD_STATUS_SUCCESS;
}

/**
 * @details     LU decomposition is applied in every diagonal matrix.
 */
ucfd_status_t lusgs_prepare(const int neles, const int nvars, const int *diag_ind, double *nnz_data)
{
    const int nvars2 = nvars * nvars;
    int idx, didx;

    // Get diagonal block and store reverse
    for (idx = 0; idx < neles; idx++)
    {
        didx = diag_ind[idx];
        ludcmp(nvars, &nnz_data[didx*nvars2]);
    }

    return UCFD_STATUS_SUCCESS;
}

/**
 * @details     This function applies preconditioner matrix into arbitrary input vector `b`.
 *              In other words, solve `Px = b`.
 */
ucfd_status_t lusgs_psolve(const int neles, const int nvars, const int *row_ptr,
                           const int *col_ind, const int *diag_ind, double *nnz_data, double *b)
{
    const int nvars2 = nvars * nvars;
    int idx, jdx, kdx, row, col;
    int dd, st, ed, cind;
    double v, arr[nvars];

    // Forward sweep : (D+L)x' = b -> x' = inv(D) * (b-Lx')
    for (idx = 0; idx < neles; idx++)
    {
        dd = diag_ind[idx];
        st = row_ptr[idx];

        // arr := b
        for (kdx = 0; kdx < nvars; kdx++)
            arr[kdx] = b[kdx + idx * nvars];

        // arr := b - Lx'
        for (jdx = st; jdx < dd; jdx++)
        {
            cind = col_ind[jdx];
            for (row = 0; row < nvars; row++)
            {
                v = 0.0;
                for (col = 0; col < nvars; col++)
                    v += nnz_data[col + row * nvars + jdx * nvars2] * b[col + cind * nvars];
                arr[row] -= v;
            }
        }

        // x' := inv(D) * (b-Lx') = inv(D) * arr
        lusubst(nvars, &nnz_data[dd*nvars2], arr);
        for (kdx=0; kdx<nvars; kdx++)
            b[kdx + idx*nvars] = arr[kdx];
    }

    // Backward sweep : (D+U)x = Dx' -> x = x' - inv(D) * Ux
    for (idx = neles - 1; idx > -1; idx--)
    {
        dd = diag_ind[idx];
        ed = row_ptr[idx + 1];

        // Initialize
        for (kdx = 0; kdx < nvars; kdx++)
            arr[kdx] = 0.0;

        // arr := Ux
        for (jdx = dd + 1; jdx < ed; jdx++)
        {
            cind = col_ind[jdx];
            for (row = 0; row < nvars; row++)
            {
                v = 0.0;
                for (col = 0; col < nvars; col++)
                    v += nnz_data[col + row * nvars + jdx * nvars2] * b[col + cind * nvars];
                arr[row] += v;
            }
        }

        // arr := inv(D) Ux
        lusubst(nvars, &nnz_data[dd*nvars2], arr);

        // b := b - inv(D) Ux
        for (kdx = 0; kdx < nvars; kdx++)
            b[kdx + idx * nvars] -= arr[kdx];
    }

    return UCFD_STATUS_SUCCESS;
}

