#ifndef PBUTILS_H
#define PBUTILS_H

void make_nei_ele(const int m, const int n, const int l, int **nei_ele);

int searcharr(int val, int *arr, int size);

void make_reordering(const int nele, const int nface, int **nei_ele, int *mapping, int *unmapping);

void make_coloring(const int m, const int n, const int l, \
                    int *icolor, int *lcolor);


#endif



