## MLAU/1D/SERIAL
Serial codes for the following one-dimensional MHD problems are available:
- `shock`... standard shock tube problems;
- `wave`... MHD wave propagation problems.

### How to run the simulation
```
>cd shock/
>make
>./a.out
```

The result is stored in `dat/`.

Simulation parameters are defined in `main.cpp`.

Users can select a numerical scheme by editting `mhd_fd4c_1d.c`:
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
  BX              FLOAT     = Array[800, 51]
  BY              DOUBLE    = Array[800, 51]
  BZ              DOUBLE    = Array[800, 51]
  DIR             STRING    = 'dat'
  DT              FLOAT     =     0.0100509
  DX              FLOAT     =    0.00124997
  EN              DOUBLE    = Array[800, 51]
  FILES           STRING    = Array[51]
  GAMMA           FLOAT     =       1.66667
  MX              DOUBLE    = Array[800, 51]
  MY              DOUBLE    = Array[800, 51]
  MZ              DOUBLE    = Array[800, 51]
  NT              LONG      =           51
  NX              LONG      =          800
  PR              DOUBLE    = Array[800, 51]
  RO              DOUBLE    = Array[800, 51]
  T               FLOAT     = Array[51]
  VX              DOUBLE    = Array[800, 51]
  VY              DOUBLE    = Array[800, 51]
  VZ              DOUBLE    = Array[800, 51]
  X               FLOAT     = Array[800]
  XOFF            INT       =        4
  ...  
```

### Read data with Python
Execute the python script `batch.py` to read the data at a particular period,
```
>python
>>>exec(open("batch.py").read())
Input data directory (Ctrl-D to exit): dat
Specity time period (0-50): 20
```
or `batch_a.py` to read the whole data.
```
>python
>>> exec(open('batch_a.py').read())
Input data directory (Ctrl-D to exit): dat
Specity time period to plot (0-50): 20
```
