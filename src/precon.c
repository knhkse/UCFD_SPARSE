#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"
#include "precon.h"


// !! Diagonal matrices are inverted !!
// !! CANNOT be applied into MKL sparse trsv function !!
ucfd_precon_status_t ilu0_prepare_bsr(const int neles, const int nvars, \
                                      const int *row_ptr, const int *col_ind, const int *diag_ind, double *nnz_data)
{
    int idx, jdx, kdx, ldx, mdx, row, col;
    int st, ed, ijdx, kjdx, ijed, kjed;
    double mat[nvars*nvars];
    int ijarr[nvars+3], kjarr[nvars+3];
    lapack_int info, *ipiv;
    
    // Allocate memory
    ipiv = (lapack_int *)malloc(sizeof(lapack_int)*nvars);

    for (idx=0; idx<neles; idx++) {
        st = row_ptr[idx];
        ed = diag_ind[idx];

        for (kdx=st; kdx<ed; kdx++) {
            ldx = diag_ind[col_ind[kdx]];
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nvars, nvars, nvars, 1.0, \
                        &nnz_data[kdx*nvars*nvars], nvars, &nnz_data[ldx*nvars*nvars], nvars, 0.0, mat, nvars);
            cblas_dcopy(nvars*nvars, mat, 1, &nnz_data[kdx*nvars*nvars], 1);

            ijdx = kdx + 1;
            kjdx = ldx + 1;
            ijed = row_ptr[idx+1];
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
                else if (col_ind[ijdx] < col_ind[kjdx]) ijdx++;
                else kjdx++;
            }

            // A[i,j] = A[i,j] - A[i,k] @ A[k,j]
            for (jdx=0; jdx<mdx; jdx++) {
                ijdx = ijarr[jdx];
                kjdx = kjarr[jdx];
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nvars, nvars, nvars, -1.0, \
                            mat, nvars, &nnz_data[kjdx*nvars*nvars], nvars, 1.0, &nnz_data[ijdx*nvars*nvars], nvars);
            }
        }
        info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, nvars, nvars, &nnz_data[ed*nvars*nvars], nvars, ipiv);
        if (info != 0) printf("Lapack dgetrf function error\n");
        info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, nvars, &nnz_data[ed*nvars*nvars], nvars, ipiv);
        if (info != 0) printf("Lapack dgetrf function error\n");
    }

    free(ipiv);

    return UCFD_PRECON_SUCCESS;
}


// * MKL sparse function ver.
ucfd_precon_status_t ilu0_sweep_bsr_mkl(sparse_matrix_t op, double *b)
{
    sparse_status_t status;
    struct matrix_descr ldescr, udescr;

    // Set matrix_description values
    // Lower-triangular matrix
    ldescr.type = SPARSE_MATRIX_TYPE_BLOCK_TRIANGULAR;
    ldescr.mode = SPARSE_FILL_MODE_LOWER;
    ldescr.diag = SPARSE_DIAG_UNIT;

    // Upper-triangular matrix
    udescr.type = SPARSE_MATRIX_TYPE_BLOCK_TRIANGULAR;
    udescr.mode = SPARSE_FILL_MODE_UPPER;
    udescr.diag = SPARSE_DIAG_NON_UNIT;

    status = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, ldescr, b, b);
    if (status != SPARSE_STATUS_SUCCESS) return UCFD_PRECON_SPARSE_FAILED;
    status = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, udescr, b, b);
    if (status != SPARSE_STATUS_SUCCESS) return UCFD_PRECON_SPARSE_FAILED;

    return UCFD_PRECON_SUCCESS;
}


// * Pure C function ver.
ucfd_precon_status_t ilu0_sweep_bsr(const int neles, const int nvars, const int *row_ptr, const int *col_ind, const int *diag_ind, \
                                    double *nnz_data, double *b)
{
    int idx, jdx, kdx, row, col;
    int dd, st, ed, cind;
    double arr[nvars];

    // Solve Lower-triangular matrix
    for (idx=0; idx<neles; idx++) {
        cblas_dcopy(nvars, &b[idx*nvars], 1, arr, 1);

        dd = diag_ind[idx];
        st = row_ptr[idx];

        for (jdx=st; jdx<dd; jdx++) {
            cind = col_ind[jdx];

            cblas_dgemv(CblasRowMajor, CblasNoTrans, nvars, nvars, -1.0, &nnz_data[jdx*nvars*nvars], nvars, \
                        &b[cind*nvars], 1, 1.0, arr, 1);
        }

        cblas_dcopy(nvars, arr, 1, &b[idx*nvars], 1);
    }

    // Solve Upper-triangular matrix
    for (idx=neles-1; idx>-1; idx--) {
        cblas_dcopy(nvars, &b[idx*nvars], 1, arr, 1);

        dd = diag_ind[idx];
        ed = row_ptr[idx+1];

        for (jdx=(dd+1); jdx<ed; jdx++) {
            cind = col_ind[jdx];

            cblas_dgemv(CblasRowMajor, CblasNoTrans, nvars, nvars, -1.0, &nnz_data[jdx*nvars*nvars], nvars, \
                        &b[cind*nvars], 1, 1.0, arr, 1);
        }

        cblas_dgemv(CblasRowMajor, CblasNoTrans, nvars, nvars, 1.0, &nnz_data[dd*nvars*nvars], nvars, \
                    arr, 1, 1.0, &b[idx*nvars], 1);
    }

    return UCFD_PRECON_SUCCESS;
}




