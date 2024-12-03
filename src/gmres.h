// Doxygen file comments


#ifndef GMRES_H
#define GMRES_H
#include "mkl.h"
#include "types.h"


ucfd_status_t serial_gmres(sparse_matrix_t op, ucfd_precon_type_t precon_type, sparse_matrix_t precon, const int n, const int m, \
                           const double tol, const double itmax, double *dub, double *rhsb, double *H, double *V, double *g, double *y, double *w, double *r);

ucfd_status_t serial_gmres2(sparse_matrix_t op, ucfd_precon_type_t precon_type, const int neles, const int nvars, const int m, \
                            const int *row_ptr, const int *col_ind, const int *diag_ind, double *pre_nnz_data, \
                            const double tol, const double itmax, double *dub, double *rhsb, double *H, double *V, double *g, double *y, double *w, double *r);


ucfd_status_t single_gmres(sparse_matrix_t op, ucfd_precon_type_t precon_type, sparse_matrix_t precon, const int n, const int m, \
                           const double tol, const double itmax, double *dub, double *rhsb, double *H, double *V, double *g, double *y, double *w, double *r);

#endif // GMRES_H