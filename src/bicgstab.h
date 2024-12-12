// Doxygen file comments

#ifndef BICGSTAB_H
#define BICGSTAB_H
#include "mkl.h"
#include "types.h"


ucfd_status_t serial_bicgstab_mkl(sparse_matrix_t op, ucfd_precon_type_t precon_type, sparse_matrix_t precon, const int n, \
                                  const double tol, const double itmax, double *dub, double *rhsb, double *r, double *p, double *v, double *s, double *t);


ucfd_status_t serial_bicgstab(sparse_matrix_t op, ucfd_precon_type_t precon_type, const int neles, const int nvars, \
                              const int *row_ptr, const int *col_ind, const int *diag_ind, double *pre_nnz_data, \
                              const double tol, const double itmax, double *dub, double *rhsb, double *r, double *p, double *v, double *s, double *t);




#endif // BICGSTAB_H