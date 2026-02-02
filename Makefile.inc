## Intel oneAPI (Recommended)
CC = icx
CFLAGS = -Wall -O3 -ffast-math -qopenmp

## GCC or Clang
# CC = mpicc
# CFLAGS = -Wall -O3 -ffast-math -fopenmp

## LU-SGS Library Variables
NVARS = 7
NFVARS = 5
NTURBVARS = 2
NDIMS = 3

## Constants (Optional)
#GAMMA = 1.4
#PMIN = 1e-3
#BETAST = 0.09

## Data types (Optional)
# FLOAT32 = 0
# INT64 = 0


## Krylov Library Variables
BLOCKSIZE = $(NVARS)


