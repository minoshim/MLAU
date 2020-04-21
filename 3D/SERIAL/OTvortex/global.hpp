#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "myfunc.h"

#define XMESH (200)		// Number of cells in X domain
#define YMESH (200)		// Number of cells in Y domain
#define ZMESH (8)		// Number of cells in Z domain
#define XOFF (4)		// Number of ghost cells in each X side.
#define YOFF (4)		// Number of ghost cells in each Y side.
#define ZOFF (4)		// Number of ghost cells in each Z side.

#define CFLCHECK (0)		// Flag to modify dt at every step

namespace global
{
  // Universal parameters (fixed)
  const int nx=XMESH+2*XOFF,ny=YMESH+2*YOFF,nz=ZMESH+2*ZOFF;
  const double pi=4.0*atan(1.0);
  const double dtor=pi/180.;

  // Parameters
  int nrec=150;			// Number of iterations for output
  int nmax=nrec*5;		// Number of maximum iteration
  const double lx=2.0*pi;	// Spatial domain in X
  const double ly=2.0*pi;	// Spatial domain in Y
  const double lz=2.0*pi;	// Spatial domain in Z
  const double dx=lx/XMESH;	// Spatial width in X
  const double dy=ly/YMESH;	// Spatial width in Y
  const double dz=lz/ZMESH;	// Spatial width in Z
  const double idx=1.0/dx;
  const double idy=1.0/dy;
  const double idz=1.0/dz;
  const double gam=5.0/3.0;	// Specific heat ratio
  const double dr=min3(dx,dy,dz);
  const double cfl=0.4;		// CFL number
  const char fildir[]="dat/";	// Directory for file output
  double vfast=3.0;		// Rough estimate of vfast+vbulk for OTvortex
  double dt=cfl*dr/vfast;	// Time step

  // Variables
  double *x,*y,*z;
  double *ro,*mx,*my,*mz,*en,*bx,*by,*bz;
}
