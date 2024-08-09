#ifndef ARRAYS_H
#define ARRAYS_H

void **malloc_2d(const size_t rows, const size_t cols, const size_t T);

void ***malloc_3d(const size_t rows, const size_t cols, const size_t depth, const size_t T);

void dealloc_2d(void **mat);

void dealloc_3d(void ***mat);



#endif