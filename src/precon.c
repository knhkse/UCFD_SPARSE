#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"
#include "precon.h"


// !! Diagonal matrices are inverted !!
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


ucfd_precon_status_t ilu0_prepare_bsr_raw(const int neles, const int nvars, \
                                          const int *row_ptr, const int *col_ind, const int *diag_ind, double *nnz_data)
{
    int idx, jdx, kdx, ldx, mdx, row, col, ele;
    int st, ed, ijdx, kjdx, ijed, kjed;
    double v, mat[nvars*nvars];
    int ijarr[nvars+3], kjarr[nvars+3];
    int nvars2 = nvars*nvars;
    lapack_int info, *ipiv;

    // Allocate memory
    ipiv = (lapack_int *)malloc(sizeof(lapack_int)*nvars);

    for (idx=0; idx<neles; idx++) {
        st = row_ptr[idx];
        ed = diag_ind[idx];

        for (kdx=st; kdx<ed; kdx++) {
            ldx = diag_ind[col_ind[kdx]];
            for (row=0; row<nvars; row++) {
                for (col=0; col<nvars; col++) {
                    v = 0.0;
                    for (ele=0; ele<nvars; ele++)
                        v += nnz_data[ele+row*nvars+kdx*nvars2] * nnz_data[col+ele*nvars+ldx*nvars2];
                    mat[col+row*nvars] = v;
                }
            }

            for (row=0; row<nvars; row++) {
                for (col=0; col<nvars; col++)
                    nnz_data[col+row*nvars+kdx*nvars2] = mat[col+row*nvars];
            }

            ijdx = kdx+1;
            kjdx = ldx+1;
            ijed = row_ptr[idx+1];
            kjed = row_ptr[col_ind[kdx]+1];
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

            // A[i,j] -= A[i,k] @ A[k,j]
            for (jdx=0; jdx<mdx; jdx++) {
                ijdx = ijarr[jdx];
                kjdx = kjarr[jdx];
                for (row=0; row<nvars; row++) {
                    for (col=0; col<nvars; col++) {
                        v = 0.0;
                        for (ele=0; ele<nvars; ele++)
                            v += mat[ele+row*nvars] * nnz_data[col+ele*nvars+kjdx*nvars2];
                        nnz_data[col+row*nvars+ijdx*nvars2] -= v;
                    }
                }
            }
        }
        // Inverse current row diagonal matrix
        info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, nvars, nvars, &nnz_data[ed*nvars*nvars], nvars, ipiv);
        if (info != 0) printf("Lapack dgetrf function error\n");
        info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, nvars, &nnz_data[ed*nvars*nvars], nvars, ipiv);
        if (info != 0) printf("Lapack dgetrf function error\n");
    }
    free(ipiv);

    return UCFD_PRECON_SUCCESS;
}


ucfd_precon_status_t ilu0_sweep_bsr(const int neles, const int nvars, const int *row_ptr, const int *col_ind, const int *diag_ind, \
                                        double *nnz_data, double *b)
{
    const int nvars2 = nvars*nvars;
    int idx, jdx, kdx, row, col;
    int dd, st, ed, cind;
    double v, arr[nvars];

    for (idx=0; idx<neles; idx++) {
        dd = diag_ind[idx];
        st = row_ptr[idx];

        // Initialize arr
        for (kdx=0; kdx<nvars; kdx++) arr[kdx] = b[kdx+idx*nvars];

        for (jdx=st; jdx<dd; jdx++) {
            cind = col_ind[jdx];

            for (row=0; row<nvars; row++) {
                v = 0.0;
                for (col=0; col<nvars; col++)
                    v += nnz_data[col+row*nvars+jdx*nvars2] * b[col+cind*nvars];
                arr[row] -= v;
            }
        }
        
        for (kdx=0; kdx<nvars; kdx++)
            b[kdx+idx*nvars] = arr[kdx];
    }

    for (idx=neles-1; idx>-1; idx--) {
        dd = diag_ind[idx];
        ed = row_ptr[idx+1];

        // Initialize
        for (kdx=0; kdx<nvars; kdx++) arr[kdx] = b[kdx+idx*nvars];

        for (jdx=dd+1; jdx<ed; jdx++) {
            cind = col_ind[jdx];

            for (row=0; row<nvars; row++) {
                v = 0.0;
                for (col=0; col<nvars; col++)
                    v += nnz_data[col+row*nvars+jdx*nvars2] * b[col+cind*nvars];
                arr[row] -= v;
            }
        }

        // Multiply diagonal inverse
        for (row=0; row<nvars; row++) {
            v = 0.0;
            for (col=0; col<nvars; col++)
                v += nnz_data[col+row*nvars+dd*nvars2] * arr[col];
            b[row+idx*nvars] = v;
        }
    }

    return UCFD_PRECON_SUCCESS;
}

