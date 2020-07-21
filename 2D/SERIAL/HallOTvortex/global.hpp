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
#define XOFF (4)		// Number of ghost cells in each X side.
#define YOFF (4)		// Number of ghost cells in each Y side.

#define CFLCHECK (1)		// Flag to modify dt at every step
#define HALL (1)		// Flag for Hall-MHD

namespace global
{
  // Universal parameters (fixed)
  const int nx=XMESH+2*XOFF,ny=YMESH+2*YOFF;
  const double pi=4.0*atan(1.0);
  const double dtor=pi/180.;

  // Parameters
  int nmax=10000;		// Number of maximum iteration
  const double dtrec=0.1;  	// Time step for output
  double trec=dtrec;
  const double tmax=1.01;	// Maximum simulation time
  
  const double lx=1.0;	// Spatial domain in X
  const double ly=1.0;	// Spatial domain in Y
  const double dx=lx/XMESH;	// Spatial width in X
  const double dy=ly/YMESH;	// Spatial width in Y
  const double idx=1.0/dx;
  const double idy=1.0/dy;
  const double gam=5.0/3.0;	// Specific heat ratio
  const double dr=((dx < dy)?dx:dy);
  const double cfl=0.2;		// CFL number
  const double ro0=gam*gam;	// Initial density
  const double pr0=gam;		// Initial pressure
  const double v0=1.0;		// Initial velocity
  const double b0=1.0;		// Initial magnetic field
  const char fildir[]="dat/";	// Directory for file output

  // Hall parameters
  const double iner_p=10.0*dr;	// Ion inertia length
  const double eta_h=HALL*iner_p*sqrt(ro0); // Hall resistivity

  double vfast=1.0;
  double vw=eta_h*pi/dr;	// Whistler phase velocity
  double dt=cfl*dr/vfast;	// Time step
  double dtw=cfl*dr/(vfast+vw);	// Time step for Hall term

  // Variables
  double *x,*y;
  double *ro,*mx,*my,*mz,*en,*bx,*by,*bz;
}
