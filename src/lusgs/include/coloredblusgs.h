#ifndef COLOREDBLUSGS_H
#define COLOREDBLUSGS_H

void ns_parallel_pre_blusgs(int neles, int nfvars, int nface, double factor, \
                            double *fnorm_vol, double *dt, double *diag, double *fjmat);


void rans_parallel_pre_blusgs(int neles, int nvars, int nfvars, int nface, double factor, double betast, \
                              double *fnorm_vol, double *uptsb, double *dt, double *tdiag, double *tjmat, double *dsrc);


void ns_parallel_block_sweep(int n0, int ne, int neles, int nfvars, int nface, \
                             int *nei_ele, int *icolor, int *lcolor, double *fnorm_vol, \
                             double *rhsb, double *dub, double *diag, double *fjmat);


void rans_parallel_block_sweep(int n0, int ne, int neles, int nvars, int nfvars, int nface, \
                             int *nei_ele, int *icolor, int *lcolor, double *fnorm_vol, \
                             double *rhsb, double *dub, double *tdiag, double *tjmat);


void parallel_update(int neles, int nvars, double *uptsb, double *dub, double *subres);



#endif