/** ======================================================================================================================
 * @file        krylov.c
 * @brief       Source file of Krylov subspace methods
 * @details     There are two Krylov subspace methods, GMRES and BiCGstab methods.  
 *                
 *              (1) GMRES (Generalized Minimal RESidual)  
 *              GMRES method is introduced by Saad and Schultz, which is the one of the Krylov subspace methods.
 *              GMRES algorithm consists of sparse matrix-vector multiplication, vector-vector product, L-2 norm,
 *              and these BLAS functions are implemented by using `Intel MKL` library.
 *              Restart version is used.  
 *                
 *              (2) BiCGstab (Bi-Conjugate Gradient stabilized)  
 *              Bi-CGSTAB method is introduced by H.A. van der Vorst, which is one of the Krylov subspace methods.
 *              It takes for q_k a product of appropriate 1-step Minimal Residual polynomials so that can converge
 *              rather smoothly than bi-conjugate gradient (Bi-CG) method.
 *              However, its computation process is slightly more expensive than CGS.
 *
 * @note        For system matrix A, any sparse matrix storage format (CSR, CSC, BSR) is available.
 *              However, Preconditioner arrays must be stored as `Block Sparse Row(BSR)` format.
 *
 * @author
 *              - Namhyoung Kim (knhkse@inha.edu), Department of Aerospace Engineering, Inha University
 *              - Jin Seok Park (jinseok.park@inha.ac.kr), Department of Aerospace Engineering, Inha University
 *
 * @date        Dec 2024
 * @version     1.0
 * @par         Copyright
 *              Copyright (c) 2024, Namhyoung Kim and Jin Seok Park, Inha University, All rights reserved.
 * @par         License
 *              This project is release under the terms of the MIT License (see LICENSE file).
 *
 * =======================================================================================================================
 */

#include <stdio.h>
#include <math.h>
#include "mkl.h"
#include "krylov.h"
#include "precon.h"

#define eps 2.22e-16 // Machine epsilon of double precision

/**
 * @details     Overall process of GMRES routine.
 *              Outer iteration ends when solution is converged or reached in maximum iteration number.
 *              UCFD_STATUS_CONVERGED is returned when L-2 norm of the residual vector becomes smaller than `tol`,
 *              and UCFD_STATUS_NOT_CONVERGED is returned when maximum iteration finished.
 */
ucfd_status_t serial_gmres(sparse_matrix_t op, ucfd_precon_type_t precon_type, const int neles, const int nvars, const int m,
                           const int *row_ptr, const int *col_ind, const int *diag_ind, double *pre_nnz_data,
                           const double tol, const double itmax, double *x, double *b, double *H, double *V, double *g, double *y, double *w, double *r)
{
    int it, i, j;
    const int n = neles * nvars;
    double beta, tmp, c, s, h1, h2, rr;
    sparse_status_t status;

    // Sparse matrix description for system matrix A(op)
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    descr.mode = 0;
    descr.diag = 0;

    // Computes residual : r := b - A @ x
    cblas_dcopy(n, b, 1, r, 1);
    status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, op, descr, x, 1.0, r);
    if (status != SPARSE_STATUS_SUCCESS)
    {
        printf("MKL Sparse matrix - vector multiplication error\n");
        return UCFD_STATUS_ERROR;
    }

    for (it = 0; it < itmax; it++)
    {
        beta = cblas_dnrm2(n, r, 1);
        if (beta < tol)
            return UCFD_STATUS_CONVERGED;
        
        // Left-preconditioning
        if (precon_type == BILU)
            bilu_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, r);
        else if (precon_type == LUSGS)
            lusgs_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, r);

        beta = cblas_dnrm2(n, r, 1);
        y[0] = beta;

        // V = r/beta
        #pragma omp parallel for
        for (i=0; i<n; i++) {
            V[i] = r[i]/beta;
        }

        for (j = 0; j < m; j++)
        {
            status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &V[j * n], 0.0, w);
            if (status != SPARSE_STATUS_SUCCESS)
            {
                printf("MKL Sparse matrix - vector multiplication error\n");
                return UCFD_STATUS_ERROR;
            }

            // TODO : Use FFI
            if (precon_type == BILU)
                bilu_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, w);
            else if (precon_type == LUSGS)
                lusgs_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, w);

            // Arnoldi iteration
            for (i = 0; i < (j + 1); i++)
            {
                tmp = cblas_ddot(n, w, 1, &V[i * n], 1);
                H[j + m * i] = tmp;
                cblas_daxpy(n, -tmp, &V[i * n], 1, w, 1);
            }

            tmp = cblas_dnrm2(n, w, 1);
            H[j + (j + 1) * m] = tmp;
            
            #pragma omp parallel for
            for (i=0; i<n; i++) {
                V[(j+1)*n+i] = w[i]/tmp;
            }

            // Givens Rotation
            for (i = 0; i < j; i++)
            {
                c = g[i * 2];     // g[i, 0]
                s = g[i * 2 + 1]; // g[i, 1]
                h1 = H[j + i * m];
                h2 = H[j + (i + 1) * m];
                H[j + i * m] = c * h1 - s * h2;
                H[j + (i + 1) * m] = s * h1 + c * h2;
            }

            h1 = H[j * (m + 1)];
            h2 = H[j + (j + 1) * m];
            rr = sqrt(h1 * h1 + h2 * h2);
            c = h1 / rr;
            s = -h2 / rr;
            H[j * (m + 1)] = c * h1 - s * h2;
            g[j * 2] = c;
            g[j * 2 + 1] = s;

            // Modify e1 vector
            y[j + 1] = y[j];
            y[j] *= c;
            y[j + 1] *= s;
        }

        // Back substitution
        cblas_dtrsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, m, H, m, y, 1);

        // Update
        cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, V, n, y, 1, 1.0, x, 1);

        // Computes next iteration residual
        cblas_dcopy(n, b, 1, r, 1);
        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, op, descr, x, 1.0, r);
        if (status != SPARSE_STATUS_SUCCESS)
            return UCFD_STATUS_ERROR;
    }

    return UCFD_STATUS_NOT_CONVERGED;
}

/**
 * @details     Single iteration of GMRES.
 * @note        Residual array `r` must be initialized, `r := b - A @ x`.
 */
ucfd_status_t step_gmres(sparse_matrix_t op, ucfd_precon_type_t precon_type, const int neles, const int nvars, const int m,
                         const int *row_ptr, const int *col_ind, const int *diag_ind, double *pre_nnz_data,
                         double *x, double *b, double *H, double *V, double *g, double *y, double *w, double *r)
{
    int i, j;
    const int n = neles * nvars;
    double beta, tmp, c, s, h1, h2, rr;
    sparse_status_t status;

    // Sparse matrix description for system matrix A(op)
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    descr.mode = 0;
    descr.diag = 0;

    if (precon_type == BILU)
        bilu_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, r);
    else if (precon_type == LUSGS)
        lusgs_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, r);
    
    beta = cblas_dnrm2(n, r, 1);
    y[0] = beta;

    // V = r/beta
    #pragma omp parallel for
    for (i=0; i<n; i++) {
        V[i] = r[i]/beta;
    }

    for (j = 0; j < m; j++)
    {
        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &V[j * n], 0.0, w);
        if (status != SPARSE_STATUS_SUCCESS)
            return UCFD_STATUS_ERROR;

        if (precon_type == BILU)
            bilu_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, w);
        else if (precon_type == LUSGS)
            lusgs_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, w);

        // Arnoldi iteration
        for (i = 0; i < (j + 1); i++)
        {
            tmp = cblas_ddot(n, w, 1, &V[i * n], 1);
            H[j + m * i] = tmp;
            cblas_daxpy(n, -tmp, &V[i * n], 1, w, 1);
        }

        tmp = cblas_dnrm2(n, w, 1);
        H[j + (j + 1) * m] = tmp;

        #pragma omp parallel for
        for (i=0; i<n; i++) {
            V[(j+1)*n+i] = w[i]/tmp;
        }

        // Givens Rotation
        for (i = 0; i < j; i++)
        {
            c = g[i * 2];     // g[i, 0]
            s = g[i * 2 + 1]; // g[i, 1]
            h1 = H[j + i * m];
            h2 = H[j + (i + 1) * m];
            H[j + i * m] = c * h1 - s * h2;
            H[j + (i + 1) * m] = s * h1 + c * h2;
        }

        h1 = H[j * (m + 1)];
        h2 = H[j + (j + 1) * m];
        rr = sqrt(h1 * h1 + h2 * h2);
        c = h1 / rr;
        s = -h2 / rr;
        H[j * (m + 1)] = c * h1 - s * h2;
        g[j * 2] = c;
        g[j * 2 + 1] = s;

        // Modify e1 vector
        y[j + 1] = y[j];
        y[j] *= c;
        y[j + 1] *= s;
    }

    // Back substitution
    cblas_dtrsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, m, H, m, y, 1);

    // Update
    cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, V, n, y, 1, 1.0, x, 1);

    // Computes next iteration residual
    cblas_dcopy(n, b, 1, r, 1);
    status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, op, descr, x, 1.0, r);
    if (status != SPARSE_STATUS_SUCCESS)
        return UCFD_STATUS_ERROR;

    return UCFD_STATUS_SUCCESS;
}

/**
 * @details     Overall process of BiCGstab routine.
 *              Outer iteration ends when solution is converged or reached in maximum iteration number.
 *              UCFD_STATUS_CONVERGED is returned when L-2 norm of the residual vector becomes smaller than `tol`,
 *              and UCFD_STATUS_NOT_CONVERGED is returned when maximum iteration finished.
 */
ucfd_status_t serial_bicgstab(sparse_matrix_t op, ucfd_precon_type_t precon_type, const int neles, const int nvars,
                              const int *row_ptr, const int *col_ind, const int *diag_ind, double *pre_nnz_data,
                              const double tol, const double itmax, double *x, double *b, double *r, double *p, double *v, double *s, double *t)
{
    int it;
    const int n = neles * nvars;
    double rho, rhoprev, alpha, beta, omega, resid;
    double rv, ts, tt;
    sparse_status_t status;

    // Sparse matrix description for system matrix A(op)
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    descr.mode = 0;
    descr.diag = 0;

    // Computes residual : r := b - A @ x
    cblas_dcopy(n, b, 1, r, 1);
    status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, op, descr, x, 1.0, r);
    if (status != SPARSE_STATUS_SUCCESS)
    {
        printf("MKL Sparse matrix - vector multiplication error\n");
        return UCFD_STATUS_ERROR;
    }

    // Choose r\tilde as r
    /* r[:n] = r
    r[n:] = r\tilde */
    cblas_dcopy(n, r, 1, &r[n], 1);

    // Outer iteration
    for (it = 0; it < itmax; it++)
    {
        resid = cblas_dnrm2(n, r, 1);
        printf("%d : %.10f\n", it, resid);
        if (resid < tol)
            return UCFD_STATUS_CONVERGED;

        rho = cblas_ddot(n, r, 1, &r[n], 1);

        // Rho breakdown
        if (fabs(rho) < eps)
            return UCFD_STATUS_RHO_BREAKDOWN;

        if (it == 0)
            cblas_dcopy(n, r, 1, p, 1);
        else
        {
            beta = (rho / rhoprev) * (alpha / omega);

            // p = r + beta*(p - omega*v)
            cblas_daxpy(n, -omega, v, 1, p, 1);
            cblas_dscal(n, beta, p, 1);
            cblas_daxpy(n, 1.0, r, 1, p, 1);
        }

        // phat = inv(M) @ p
        cblas_dcopy(n, p, 1, &p[n], 1);

        // Apply preconditioner to p[n]~
        if (precon_type == BILU)
            bilu_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, &p[n]);
        else if (precon_type == LUSGS)
            lusgs_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, &p[n]);

        // v = A @ phat
        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &p[0], 0.0, v);
        if (status != SPARSE_STATUS_SUCCESS)
        {
            printf("MKL Sparse matrix - vector multiplication error\n");
            return UCFD_STATUS_ERROR;
        }

        rv = cblas_ddot(n, v, 1, &r[n], 1);
        alpha = rho / rv;

        // s = r - alpha*v
        // 1) r := r - alpha*v (r = s)
        cblas_daxpy(n, -alpha, v, 1, r, 1);
        // 2) s := r
        cblas_dcopy(n, r, 1, s, 1);

        // x := x + alpha*phat
        cblas_daxpy(n, alpha, &p[n], 1, x, 1);
        resid = cblas_dnrm2(n, s, 1);
        if (resid < tol)
            return UCFD_STATUS_CONVERGED;

        // shat = inv(M) @ s
        cblas_dcopy(n, s, 1, &s[n], 1);

        if (precon_type == BILU)
            bilu_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, &s[n]);
        else if (precon_type == LUSGS)
            lusgs_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, &s[n]);

        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &s[0], 0.0, t);
        if (status != SPARSE_STATUS_SUCCESS)
        {
            printf("MKL Sparse matrix - vector multiplication error\n");
            return UCFD_STATUS_ERROR;
        }

        ts = cblas_ddot(n, t, 1, s, 1);
        tt = cblas_ddot(n, t, 1, t, 1);
        omega = ts / tt;

        // Update solution
        cblas_daxpy(n, omega, &s[0], 1, x, 1);

        // r = s - omega*t
        cblas_daxpy(n, -omega, t, 1, r, 1);

        // Update rho
        rhoprev = rho;
    }

    return UCFD_STATUS_ERROR;
}
