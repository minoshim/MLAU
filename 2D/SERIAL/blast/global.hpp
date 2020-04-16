#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "myfunc.h"

#define XMESH (1024)		// Number of cells in X domain
#define YMESH (1024)		// Number of cells in Y domain
#define XOFF (4)		// Number of ghost cells in each X side.
#define YOFF (4)		// Number of ghost cells in each Y side.

#define CFLCHECK (1)		// Flag to modify dt at every step

namespace global
{
  // Universal parameters (fixed)
  const int nx=XMESH+2*XOFF,ny=YMESH+2*YOFF;
  const double pi=4.0*atan(1.0);
  const double dtor=pi/180.;

  // Parameters
  int nmax=2000;		// Number of maximum iteration
  const double dtrec=0.01;	// Time step for output  
  double trec=dtrec;

  const double lx=4.0;	// Spatial domain in X
  const double ly=lx*(double)YMESH/XMESH;// Spatial domain in Y
  const double dx=lx/XMESH;	// Spatial width in X
  const double dy=ly/YMESH;	// Spatial width in Y
  const double idx=1.0/dx;
  const double idy=1.0/dy;
  const double gam=5.0/3.0;	// Specific heat ratio
  const double dr=((dx < dy)?dx:dy);
  const double cfl=0.4;		// CFL number
  const char fildir[]="dat/";	// Directory for file output

  const double r_0=0.125;	// Radius of imposed high-P cylinder
  const double ro0=1e0;		// Ambient density
  const double ro1=1e0;		// Density in cylinder
  const double pr0=1e0;		// Ambient pressure
  const double pr1=1e2;		// Pressure in cylinder
  const double b_0=10.0;	// Ambient |B|
  const double ban=30.0;	// B angle relative to x axis 

  double vfast=1.0;
  double dt=cfl*dr/vfast;	// Time step

  // Variables
  double *x,*y;
  double *ro,*mx,*my,*mz,*en,*bx,*by,*bz;
}
