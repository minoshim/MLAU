## MLAU/2D/MPI
MPI parallel codes for the following two-dimensional MHD problems are available:
- `KHI` ... Kelvin-Helmholtz instability;
- `MRX` ... Magnetic reconnection;
- `OTvortex` ... Orszag-Tang vortex problem;
- `RMI` ... Richtmyer-Meshkov instability.

### How to run the simulation
```
>cd OTvortex/
>make
>mpiexec -np 16 -genv OMP_NUM_THREADS 2 ./a.out #for MPICH users
>mpiexec -np 16 -x OMP_NUM_THREADS=2 ./a.out    #for OpenMPI users
```

Here `16` is the number of MPI processes and `2` is the number of OpenMP threads, thus 32 CPU cores are used for this example run.

The number of MPI processes should be equal to the value of `MNP` defined in `global.hpp` (otherwise, the simulation does not run).

Users can abort the run by Ctrl+C, and restart it by the same `a.out` and command.

The result is stored in `dat/`.

Simulation parameters are defined in `global.hpp`:<br>
`XMESH` and `YMESH` define the number of calculation cells in the whole X and Y domain, and they should be an integer multiple of the number of MPI processes in X and Y directions, `MNP_X` and `MNP_Y`.

Users can select a numerical scheme by editting `mhd_fd_ct_2d.c`:
- `ODR` (1, 2, 3, or 4) for spatial order;
- `R_K` (1, 2, or 3) for temporal order;
- `func_flux` (MLAU, HLLD, ROE, or LHLLD) for flux function.

### Read data with IDL

Execute the scripts `rddt.pro` and `read_mhd.pro` found in `idl/`.
```
>idl
  IDL> .r rddt
  % Compiled module: $MAIN$.
  Set data directory: dat
  % Compiled module: FILE_READS.
  % $MAIN$: Independent varaibles are loaded.
  % $MAIN$: To load MHD data, call READ_MHD.
  IDL> .r read_mhd
  % Compiled module: $MAIN$.
  Set time step (0 - 5) : 5
  Read MHD data at 5 step
  % Compiled module: BINARY_READ.
  IDL> help
  ...
  B2              FLOAT     = Array[200, 200]
  BX              FLOAT     = Array[200, 200]
  BX_CT           FLOAT     = Array[200, 200]
  BY              FLOAT     = Array[200, 200]
  BY_CT           FLOAT     = Array[200, 200]
  BZ              FLOAT     = Array[200, 200]
  DIR             STRING    = 'dat'
  ...
```

### Read data with Python

Users firtly merge the raw simulation data that is MPI-decomposed,
```
>./merge.out dat/ dat/
```

Subsequently, execute the python script `batch.py` to read the data at a particular period.
```
>python
>>>exec(open("batch.py").read())
Input data directory (Ctrl-D to exit): dat
Specity time period (0-5): 5
```
