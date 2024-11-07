# UCFD_SPARSE: Unstructured Grid Based CFD Applicable Asymmetric Sparse Matrix Numerical Library
[![](https://img.shields.io/badge/documentation-32CD32)](https://xccels.github.io/UCFD_SPARSE/)
![](https://img.shields.io/badge/license-MIT_License-yellow.svg)

UCFD_SPARSE provides LU-SGS (Lower Upper Symmetric Gauss-Seidel) type numerical methods for unstructured grid based Computational Fluid Dynamics (CFD) simulation. Since time accuracy is less important in steady flows problem, many researchers has intorduced robust and fast-converging time marching schemes. LU-SGS method is one of the widely-used method which was proposed by *[Jameson & Yoon (1987)](https://arc.aiaa.org/doi/abs/10.2514/3.9724)* to obtain steady solution, and has been successfully applied to the unstructured grid problem (*[Sharov & Nakahashi (1997)](https://arc.aiaa.org/doi/10.2514/6.1997-2102)*).  
2-dimensional or 3-dimensional problem with Euler/Navier-Stokes/RANS equations is applicable. For RANS equations, only 1-equation and 2-equations model are supported and are computed separately with Navier-Stokes equations.


# Features
Currently, UCFD_SPARSE module supports CPU processor only.

## Colored LU-SGS method
Original LU-SGS method has `data dependency`, so that it cannot be parallelized directly. Therefore, a lot of schemes have been developed to apply the LU-SGS method to shared memory parallelism. UCFD_SPARSE provides a parallelizable LU-SGS method using a `Multi-Coloring Algorithm`. The coloring algorithm divides grid cells into clusters which are independent of each other. Since every cell in a specific cluster does not have any data dependency, its material properties can be updated simultaneously.

## Array parameters
Currently, functions in LU-SGS source file is fit for [pyBaram](https://gitlab.com/aadl_inha/pyBaram) formulation. There are some array parameters to compute solution for next time step.  

- 1) `fnorm_vol` : face normalized volume array [nface, neles]  
    Functionally it is the same as the Jacobian term in a structured grid. This array can be made by multiplying the reciprocal of the cell volume by the surface magnitude.

- 2) `dt` : time step array [neles]  
    If `local time stepping` method is used, each value of the array contains the time step of the corresponding cell. Otherwise, all values of the array have the same value.

- 3) `fspr` : face spectral radius array [nface, neles]  
    The spectral radius of each cell face. A general implicit scheme uses Jacobian matrices, by contrast.

- 4) `nei_ele` : neighbor element array [nface, neles]  
    Neighbor cell index of each cell. If there is no neighbor cell at a face (boundary), the current cell index is implemented.

- 5) `vec_fnorm` : face normal vector [nface, ndims, neles]  
    Normalized vector of each cell face. The vector is always outward-directed with reference to the current cell.

More informations about function parameters are available in the [documents](https://xccels.github.io/UCFD_SPARSE/).



## Computation process
LU-SGS method consists of two main sweeps (lower and upper sweep), and each iteration step is executed with 4 main process.

- 1) LU-SGS preparation  
    In the preparation phase, diagonal matrix of the implicit operator is constructed. The diagonal matrix is slightly different from the off-diagonal matrix, including time step and cell volume. Generally, inviscid flux Jacobian term is approximated with a diffusive scheme like Rusanov flux, so that the diagonal matrix is reformulated by an identity matrix with a scale factor. This replaces the matrix inversion process with simple scalar division, which reduces computation time dramatically.

- 2) Lower sweep  
    After constructing the diagonal matrix, lower and upper sweeps are executed. The lower sweep starts from the least numbered cell to the largest numbered cell. The lower sweep results in an intermediate solution $\Delta \vec{Q}^*$, which is stored in the `dub` array.

- 3) Upper sweep  
    In contrast to the lower sweep, the upper sweep explores the cell reversely, from the largest numbered cell to the least numbered cell. The colored LU-SGS also has an almost identical procedure, reversing the order of each cell cluster to compute. The upper sweep results in the difference between the current and next time step solution $\Delta \vec{Q}$, which is stored in the `rhsb` array.

- 4) Update  
    $\Delta \vec{Q}$ is stored in the `rhsb` array because right hand side of the equations is no longer needed after the upper sweep. The solution of the next iteration step can be obtained by adding the difference to the current solution.


# Authors
- Namhyoung Kim (knhkse@inha.edu), Department of Aerospace Engineering, Inha University
- Jin Seok Park (jinseok.park@inha.ac.kr), Department of Aerospace Engineering, Inha University


# Usage

## Download
```
git clone https://github.com/xccels/UCFD_SPARSE.git
```

## Compile
- GCC or Intel oneAPI compiler for C
- Python interpreter (optional)

### Compile and build
- Build UCFD_SPARSE module
    ```
    make lib
    ```

- Build example problem after build UCFD_SPARSE
    ```
    make example
    ```

- Build all
    ```
    make all
    ```

### Compiler option
`Makefile.inc` file in the root directory defines compiler and optimization options. This can be modified depending on each user's computer setting.

### Running the example
After building the example, two executable files will be located in the `run` directory. The `.x` file can be executed as follows:

```
mpirun -np 8 ./mpi3d.x ./input.dat
```

or

```
mpirun -np 8 ./omp3d.x ./input.dat
```

Note that the total number of processors to use for MPI communication must match the total number of elements defined in the `input.dat` file. The total number of elements is computed by multiplying `npx`, `npy`, and `npz`.  
Python examples are also available. In this case, a dynamic library (.so) is used. Since the Python file should be executed in the `example` directory, absolute paths for the dynamic library and input file are required. Also, for more stability, it is recommended to pass the number of threads to the command arguments for multi-threading case. Python example can be executed as follows:

```
mpirun -np 8 python omp3d.x $(UCFD_PATH)/run/input.dat $(UCFD_PATH)/lib/liblusgs.so nthreads=8
```


# Folder structure
- `example` : Simple examples for using UCFD_SPARSE module. There are two different examples, using MPI and MPI+OpenMP hybrid case.
- `src` : Source files of UCFD_SPARSE
- `lib` : Static libraries of LU-SGS and Colored LU-SGS are created after building.
- `obj` : Object files for each source file are created after building.
- `run` : Executable file for the example problem is created after building.
- `utils` : Contains useful functions for example problem.
- `include` : Contains `const.h` header file which defines constants for phyisical properties. These values are not modified unless library is built again.


# References
[XCCELS](https://xccels.github.io/main/)
[Aerodynamic Analysis & Design Laboratory](http://aadl.inha.ac.kr)
[pyBaram](https://gitlab.com/aadl_inha/pyBaram)

