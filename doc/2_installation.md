Installation                        {#install_page}
============

[TOC]
# Download

## 1. Git clone

```
git clone https://github.com/xccels/UCFD_SPARSE.git
```

## 2. Download ZIP
Code - Download source code - zip
\image html install.png width=80%

# Compile
## Prerequisites

Prerequisites to compile UCFD_SPARSE are as follows:
- MPI (OpenMPI, INtelMPI, MPICH, etc)
- C Compiler (GCC, Intel C compiler, etc)
- Intel MKL
    - 


## Compile and build

- Designate module path
```
export UCFD_PATH=$(pwd)
```

### Build UCFD_SPARSE 
- make all library
```
make lib
```

- Build static library (lusgs.a)
```
make static
```

- Build dynamic library (lusgs.so)
```
make dynamic
```

### Build example
- Build example files after build UCFD_SPARSE
```
make example
```

- Build all
```
make all
```

- Clean build
```
make clean
```

## Configuration
Compilation configuration is needed when compile library.

### LU-SGS type library

* Data types
    - `INT64` : Set integer type as 64-bit. Default is int32.
    - `FLOTE32` : Set real type as 32-bit (single precision). Default is float64.

* Flow variables
    - `NVARS` : Number of total flow variables
    - `NFVARS` : Number of convective flow variables
    - `NTURBVARS` : Number of turbulent flow variables
    - `NDIMS` : Number of grid dimensions. Only 1, 2, or 3 is supported.
* Flow constants
    - `GAMMA` : Specific heat ratio
    - `PMIN` : Minimum pressure value
    - `BETAST` : Constant for kw-SST RANS model

### Krylov subspace method type library

* Intel MKL dependent configuration
    - `MKLROOT` : Intel MKL root directory path
    - `MKLFLAGS` : command line for implementing Intel MKL. Refer to the [Link Line Advisor](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html) for more details.
    - `BLOCK` : Block size fo

