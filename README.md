# MLAU/LHLLD
This C++ package is aimed to numerically solve compressible MHD equations with the low-dissipation quasi all-speed Riemann solvers:
- MLAU - A Multistate Low-dissipation Advection Upstream Splitting Method for Ideal Magnetohydrodynamics[^1];
- LHLLD- A low-dissipation HLLD approximate Riemann solver for a very wide range of Mach numbers[^2].

The package also adopts the following techniques for robust and accurate numerical simulations:
- up to 4th order accuracy in space and 3rd order accuracy in time;
- preservation of the solenoidal condition of the magnetic field by a well-designed Contrained Transport method[^3].

The current version supports 1D-3D space in Cartesian coordinates and hybrid MPI/OpenMP parallel computation.
 
Copyright 2020,2021 Takashi Minoshima

Contact at Takashi Minoshima <minoshim@jamstec.go.jp>

## System requirements
Following packages are required to be installed on the system:
- Git to install and update the code;
- C++ compiler (GNU, Intel);
- MPI library (MPICH, OpenMPI) to use the MPI parallel code;

To read and visualize the data, Interactive Data Language (IDL, commercial) or Python 3.X with NumPy and matplotlib are required (the latter can be downloaded from [Anaconda](https://www.anaconda.com/products/distribution))

The code is tested on Linux OSs (Ubuntu, Linux Mint, CentOS, including Windows Subsystem for Linux).

## Installation
1. Download the code from GitHub via `>git clone hppts://github.com/minoshim/MLAU`.
2. Move to the main directory `MLAU/`.
3. Check `Makefile.inc` and edit environment variables `CC`, `MPICC`, and `CFLAGS` to meet users environment.
4. Move to the problem directories (e.g., 1D/SERIAL/) and see README how to run the simulation.

Since the code is updated without notice, you may update the code via `>git pull origin main`.

## Journal sites:
- https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee (MLAU)
- https://www.sciencedirect.com/science/article/pii/S0021999121005349 (LHLLD)

## arXiv sites:
- http://arxiv.org/abs/2004.08012 (MLAU)
- https://arxiv.org/abs/2108.04991 (LHLLD)

[^1]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
[^2]: [Minoshima T. and Miyoshi T. 2021, JCP](https://www.sciencedirect.com/science/article/pii/S0021999121005349)
[^3]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)
