# MLAU/LHLLD
This C++ package is aimed to solve compressible MHD equations with the low-dissipation quasi all-speed Riemann solvers:
- MLAU - A Multistate Low-dissipation Advection Upstream Splitting Method for Ideal Magnetohydrodynamics[^1];
- LHLLD- A low-dissipation HLLD approximate Riemann solver for a very wide range of Mach numbers[^2].

The package also adopts the following techniques for robust and accurate numerical simulations:
- up to 4th order accuracy in space and 3rd order accuracy in time;
- preservation of the solenoidal condition of the magnetic field by a well-designed Contrained Transport method[^3].

The current version supports 1D-3D space in Cartesian coordinates and hybrid MPI/OpenMP parallel computation.
 
Copyright 2020,2021 Takashi Minoshima

Contact at Takashi Minoshima <minoshim@jamstec.go.jp>

## Github site:
- https://github.com/minoshim/MLAU

## Journal sites:
- https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee (MLAU)
- https://www.sciencedirect.com/science/article/pii/S0021999121005349 (LHLLD)

## arXiv sites:
- http://arxiv.org/abs/2004.08012 (MLAU)
- https://arxiv.org/abs/2108.04991 (LHLLD)

[^1]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
[^2]: [Minoshima T. and Miyoshi T. 2021, JCP](https://www.sciencedirect.com/science/article/pii/S0021999121005349)
[^3]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)
