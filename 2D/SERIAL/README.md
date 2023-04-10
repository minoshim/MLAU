## qasMHD/2D/SERIAL
Serial codes for the following two-dimensional MHD problems are available:
- `OTvortex` ... Orszag-Tang vortex problem;
- `blast` ... blast wave propagation problem;
- `loop` ... field loop advection problem.

### How to run the simulation
```
>cd OTvortex/
>make
>./a.out
```

The result is stored in `dat/`.

Simulation parameters are defined in `global.hpp`.

Users can select a numerical scheme by editting `mhd_fd_ct_2d.c`:
- `ODR` (1, 2, 3, or 4) for spatial order;
- `R_K` (1, 2, or 3) for temporal order;
- `func_flux` (MLAU, HLLD, ROE, or LHLLD) for flux function.

### Read data with IDL
Execute the script `rddt.pro` found in `idl/`.

```
>idl
  IDL> .r rddt
  % Compiled module: $MAIN$.
  Set data directory: dat
  % Compiled module: FILE_READS.
  % Compiled module: BINARY_READ.
  IDL> help
  ...
  BX              DOUBLE    = Array[200, 200, 6]
  BY              DOUBLE    = Array[200, 200, 6]
  BZ              DOUBLE    = Array[200, 200, 6]
  DIR             STRING    = 'dat'
  DT              FLOAT     =      0.628319
  DX              FLOAT     =     0.0314159
  DY              FLOAT     =     0.0314159
  EN              DOUBLE    = Array[200, 200, 6]
  FILES           STRING    = Array[6]
  GAMMA           FLOAT     =       1.66667
  MX              DOUBLE    = Array[200, 200, 6]
  MY              DOUBLE    = Array[200, 200, 6]
  MZ              DOUBLE    = Array[200, 200, 6]
  NT              LONG      =            6
  NX              LONG      =          200
  NY              LONG      =          200
  PR              DOUBLE    = Array[200, 200, 6]
  RO              DOUBLE    = Array[200, 200, 6]
  T               FLOAT     = Array[6]
  VX              DOUBLE    = Array[200, 200, 6]
  VY              DOUBLE    = Array[200, 200, 6]
  VZ              DOUBLE    = Array[200, 200, 6]
  X               FLOAT     = Array[200]
  XOFF            INT       =        4
  Y               FLOAT     = Array[200]
  YOFF            INT       =        4
```

### Read data with Python
Execute the python script `batch.py`.

```
>python
>>>exec(open("batch.py").read())
```
