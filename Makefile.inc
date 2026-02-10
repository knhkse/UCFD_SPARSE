# Common files
SRCDIR = $(UCFD_PATH)/src
CONFIG = $(SRCDIR)/config.h
TYPES = $(SRCDIR)/ucfd_types.h
IFLAGS = -I$(SRCDIR)


# Compiler setting
## Intel oneAPI (Recommended)
CC = icx
CFLAGS = -Wall -O3 -ffast-math -qopenmp

## GCC or Clang
# CC = mpicc
# CFLAGS = -Wall -O3 -ffast-math -fopenmp


# LU-SGS macros
## Equation parameters
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


# Krylov macros
## Block size
BLOCK = 4

## SDL command line
MKLFLAGS = -L${MKLROOT}/lib -lmkl_rt -lpthread -lm -ldl

