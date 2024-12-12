#ifndef PRECON_H
#define PRECON_H
#include "types.h"

ucfd_precon_status_t ilu0_prepare_bsr(const int neles, const int nvars, \
                                      const int *row_ptr, const int *col_ind, const int *diag_ind, double *nnz_data);

ucfd_precon_status_t ilu0_prepare_bsr_raw(const int neles, const int nvars, \
                                          const int *row_ptr, const int *col_ind, const int *diag_ind, double *nnz_data);

ucfd_precon_status_t ilu0_prepare_bsr_mkl(const int neles, const int nvars, \
                                          const int *row_ptr, const int *col_ind, const int *diag_ind, double *nnz_data);                                          

ucfd_precon_status_t ilu0_sweep_bsr_mkl(sparse_matrix_t op, double *b);

ucfd_precon_status_t ilu0_sweep_bsr(const int neles, const int nvars, const int *row_ptr, const int *col_ind, const int *diag_ind, \
                                    double *nnz_data, double *b);

ucfd_precon_status_t ilu0_sweep_bsr_raw(const int neles, const int nvars, const int *row_ptr, const int *col_ind, const int *diag_ind, \
                                        double *nnz_data, double *b);




#endif // PRECON.H