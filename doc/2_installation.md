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
\image html install.png width=50%

# Compile
## Prerequisites

Prerequisites to compile UCFD_SPARSE are as follows:
- MPI (OpenMPI, INtelMPI, MPICH, etc)
- C Compiler (GCC, Intel C compiler, etc)


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


