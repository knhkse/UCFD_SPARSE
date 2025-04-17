/** ======================================================================================================================
 * @file        bicgstab.c
 * @brief       Bi-Conjugate Gradient stabilized method
 * @details     Bi-CGSTAB method is introduced by H.A. van der Vorst, which is one of the Krylov subspace methods.
 *              It takes for q_k a product of appropriate 1-step Minimal Residual polynomials so that can converge
 *              rather smoothly than bi-conjugate gradient (Bi-CG) method.
 *              However, its computation process is slightly more expensive than CGS.
 * 
 * @note        Sparse matrix must be `Block Sparse Row(BSR)` format.
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
#include "bicgstab.h"
#include "precon.h"

#define eps 2.22e-16    // Machine epsilon of double precision

/**
 * @details     Overall process of BiCGstab routine.
 *              Outer iteration ends when solution is converged or reached in maximum iteration number.
 *              UCFD_STATUS_CONVERGED is returned when L-2 norm of the residual vector becomes smaller than `tol`,
 *              and UCFD_STATUS_NOT_CONVERGED is returned when maximum iteration finished.
 */
ucfd_status_t serial_bicgstab(sparse_matrix_t op, ucfd_precon_type_t precon_type, const int neles, const int nvars, \
                              const int *row_ptr, const int *col_ind, const int *diag_ind, double *pre_nnz_data, \
                              const double tol, const double itmax, double *x, double *b, double *r, double *p, double *v, double *s, double *t)
{
    int it, i, j, err;
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
    if (status != SPARSE_STATUS_SUCCESS) {
        printf("MKL Sparse matrix - vector multiplication error\n");
        return UCFD_STATUS_ERROR;
    }

    // Choose r\tilde as r
    /* r[:n] = r
       r[n:] = r\tilde */
    cblas_dcopy(n, r, 1, &r[n], 1);

    // Outer iteration
    for (it=0; it<itmax; it++) {
        resid = cblas_dnrm2(n, r, 1);
        printf("%d : %.10f\n", it, resid);
        if (resid < tol) return UCFD_STATUS_CONVERGED;

        rho = cblas_ddot(n, r, 1, &r[n], 1);
        
        // Rho breakdown
        if (fabs(rho) < eps) return UCFD_STATUS_RHO_BREAKDOWN;

        if (it == 0) cblas_dcopy(n, r, 1, p, 1);
        else {
            beta = (rho/rhoprev) * (alpha/omega);

            // p = r + beta*(p - omega*v)
            cblas_daxpy(n, -omega, v, 1, p, 1);
            cblas_dscal(n, beta, p, 1);
            cblas_daxpy(n, 1.0, r, 1, p, 1);
        }

        //phat = inv(M) @ p
        cblas_dcopy(n, p, 1, &p[n], 1);

        // Apply preconditioner to p[n]~
        if (precon_type == BILU) {
            bilu_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, &p[n]);
        }
        else if (precon_type == LUSGS) {
            lusgs_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, &p[n]);
        }

        // v = A @ phat
        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &p[0], 0.0, v);
        if (status != SPARSE_STATUS_SUCCESS) {
            printf("MKL Sparse matrix - vector multiplication error\n");
            return UCFD_STATUS_ERROR;
        }

        rv = cblas_ddot(n, v, 1, &r[n], 1);
        alpha = rho/rv;

        // s = r - alpha*v
        // 1) r := r - alpha*v (r = s)
        cblas_daxpy(n, -alpha, v, 1, r, 1);
        // 2) s := r
        cblas_dcopy(n, r, 1, s, 1);

        // x := x + alpha*phat
        cblas_daxpy(n, alpha, &p[n], 1, x, 1);
        resid = cblas_dnrm2(n, s, 1);
        if (resid < tol) return UCFD_STATUS_CONVERGED;

        // shat = inv(M) @ s
        cblas_dcopy(n, s, 1, &s[n], 1);

        if (precon_type == BILU) {
            bilu_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, &s[n]);
        }
        else if (precon_type == LUSGS) {
            lusgs_psolve(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, &s[n]);
        }

        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &s[0], 0.0, t);
        if (status != SPARSE_STATUS_SUCCESS) {
            printf("MKL Sparse matrix - vector multiplication error\n");
            return UCFD_STATUS_ERROR;
        }

        ts = cblas_ddot(n, t, 1, s, 1);
        tt = cblas_ddot(n, t, 1, t, 1);
        omega = ts/tt;

        // Update solution
        cblas_daxpy(n, omega, &s[0], 1, x, 1);

        // r = s - omega*t
        cblas_daxpy(n, -omega, t, 1, r, 1);

        // Update rho
        rhoprev = rho;
    }

    return UCFD_STATUS_ERROR;
}
