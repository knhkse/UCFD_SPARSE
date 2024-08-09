Introduction                       {#intro_page}
============

[TOC]

# Overview

Unstructured Grid Based CFD Applicable Asymmetric Sparse Matrix Numerical Library  

UCFD_SPARSE provides LU-SGS (Lower Upper Symmetric Gauss-Seidel) type numerical methods for unstructured grid based Computational Fluid Dynamics (CFD) simulation. Since time accuracy is less important in steady flows problems, many researchers have introduced robust and fast-converging time marching schemes. The LU-SGS method is one of the widely-used methods proposed by *[Jameson & Yoon (1987)](https://arc.aiaa.org/doi/abs/10.2514/3.9724)* to obtain a steady solution, and has been successfully applied to the unstructured grid problem (*[Sharov & Nakahashi (1997)](https://arc.aiaa.org/doi/10.2514/6.1997-2102)*).  

This library can be imported as a time integration process into CFD simulation code based on Finite Volume Method. Steady flow problems with Euler/Navier-Stokes/RANS equations are applicable, and both 2-dimensional and 3-dimensional grid are supported.  

The main process called `time integration` using LU-SGS method consists of the following 4 steps.  

1) Preparation : {serial/parallel}\_pre\_lusgs  
	In preparation phase, diagonal matrix of implicit operator is constructed. Diagonal matrix is slightly different with off-diagonal matrix, including time step and cell volume. Generally, inviscid flux Jacobian term is approximated with diffusive scheme like Rusanov flux, so that the diagonal matrix is reformulated by identity matrix with a scale factor. This replaces the matrix inversion process with simple scalar division, which reduces computation time dramatically.  

2) Lower Sweep : {ns/rans}\_{serial/parallel}\_lower\_sweep  
	After constructing diagonal matrix, lower and upper sweep are executed. Since lower sweep starts from least numbered cell to largest numbered cell, it is also called as `Forward Sweep`. Lower sweep results in the intermediate solution \f$\Delta \vec{Q}^*\f$, which is stored in `dub` array. In this process, only lower neighbor cells already updated are required.  

3) Upper Sweep : {ns/rans}\_{serial/parallel}\_upper\_sweep  
	In contrast with the lower sweep, upper sweep explores cell reversely, from largest numbered cell to least numbered cell so that it is called as `Backward Sweep`. Colored LU-SGS also has the almost same procedure, reversing order of each cell cluster to compute. Upper sweep results in the difference between current and next time step solution \f$\Delta \vec{Q}\f$. In this process, values in `rhsb` array are overwritten by the \f$\Delta \vec{Q}\f$ because right-hand-side array is no longer needed when upper sweep executes. Similar to the lower sweep step, only upper neighbor cells already updated are required.  

4) Update : {serial/parallel}\_update  
	Current solution array is updated to the next time step solution. Solution of the next iteration step can be obtained by adding difference with the current solution. After updating, arrays to compute next time step solution (right-hand-side, time step, face spectral radius, etc) and total residual will be re-calculated. Process between update and preparation step depends on the types of which CFD solver is used, but generally second-order-accurate numerical method is applied into space discretization to capture shock wave well or get more precise flow results.  

Generally, iterative time-stepping schemes like Block Jacobi, Gauss-Seidel method need `sub-iteration` steps to make solution difference array converge in current time step. However, sub-iteration steps can be removed when using LU-SGS method due to its remarkable robustness and stability.

# Features

## LU-SGS method
\image html lusgs.png width=50%

- There are two options in UCFD_SPARSE module.
	- {ns/rans} : Navier-Stokes equations / RANS equations
	- {serial/parallel} : LU-SGS / Colored LU-SGS
		- If OpenMP is enabled, both serial and parallel options are available.

- {ns/rans}\_{serial/parallel}\_lower\_sweep : A function to execute lower sweep
- {ns/rans}\_{serial/parallel}\_upper\_sweep : A function to execute upper sweep

## RANS equations
RANS (Reynolds-Averaged Navier-Stokes) equations are supported by using integral formulation.

$$
\begin{align}
\frac{\partial}{\partial t} \int_V {\vec{W}_T} \; dV + \oint_{\partial V} {(\vec{F}_{c_T} - \vec{F}_{v_T})} \; dS = \int_V {\vec{\mathbb{S}}_T \; dV}
\end{align}
$$

$$
\begin{align}
\vec{W}_T & : \text{Conservative variables} \\\
\vec{F}_{c_T} & : \text{Convective flux} \\\
\vec{F}_{v_T} & : \text{Viscous flux} \\\
\vec{\mathbb{S}}_T & : \text{Source term} \\\
\end{align}
$$

- Subscription \f$T\f$ means `Turbulent`.  
- Equation above is almost identical with the Euler or Navier-Stokes equations.
- flux.c : Only uses `convective flux` term based on `Rusanov flux scheme` *[Rusanov, V.V. (1962)](https://www.sciencedirect.com/science/article/abs/pii/0041555362900629)*.
- Convective flux of the RANS equations is calculated by multiplying contravariant velocity and turbulent conservative variables.


## Shared Memory Parallelism

- Colored LU-SGS method supports SMP(OpenMP) computing.
- Coloring array must be constructed before time integration begins.
- Refer to the `coloredlusgs.c` file in `src` folder.


# Authors
- Namhyoung Kim (knhkse@inha.edu), Department of Aerospace Engineering, Inha University
- Jin Seok Park (jinseok.park@inha.ac.kr), Department of Aerospace Engineering, Inha University


# Versions
- 1.0
	- Initial distribution