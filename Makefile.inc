# Build directory
BUILD_DIR = $(UCFDPATH)/build

# Optimization flags
OPTFLAGS = -O3 -ffast-math

# Multi-threading
USE_OMP = 1

# Use Intel MKL (Recommended)
USE_MKL = 1

# Compiler setting
## Intel oneAPI (Recommended)
CC = mpiicx
CXX = mpiicpx

## GCC or Clang
# CC = mpicc
# CXX = mpicxx

## NVCC
NVCC = nvcc
# NVCC = /usr/local/cuda/bin/nvcc


## Constants (Optional)
# GAMMA = 1.4
# PMIN = 1e-3
# BETAST = 0.09

## Data types (Optional)
# USE_FLOAT32 = 0
# USE_INT64 = 0
