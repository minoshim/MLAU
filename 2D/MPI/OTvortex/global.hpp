#include "mpi.h"
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

#define MNP_X (4)		// Number of MPI processes in X
#define MNP_Y (4)		// Number of MPI processes in Y
#define MNP (MNP_X*MNP_Y)	// Total number of MPI processes

#define CFLCHECK (0)		// Flag to modify dt at every step

int dnxs[8]={0,0,0,0,0,0,0,0};
int dnys[8]={0,0,0,0,0,0,0,0};
// Boundary condition flag for ro,mx,my,mz,en,bx,by,bz (be sure of variable order)
// 0=Periodic, +1=Neumann, -1=Dirichlet
int stxs[8]={0,0,0,0,0,1,0,0};
int stys[8]={0,0,0,0,0,0,1,0};
// Staggered grid flag for ro,mx,my,mz,en,bx,by,bz (be sure of variable order)
// Do NOT change

namespace global
{
  // Universal parameters (fixed)
  const int xoff=4,yoff=xoff;	// Number of ghost cells in each side
  const int nx=XMESH/MNP_X+2*xoff,ny=YMESH/MNP_Y+2*yoff;
  const int mpi_numx=MNP_X;
  const int mpi_numy=MNP_Y;
  const double pi=4.0*atan(1.0);
  const double dtor=pi/180.;

  // Parameters
  int nrec=150;			// Number of steps for output
  int nmax=nrec*5;		// Number of maximum iteration
  const double dtrec=0.2*pi; 	// Time step for output  
  double trec=dtrec;

  const double lx=2.0*pi;	// Spatial domain in X
  const double ly=2.0*pi;	// Spatial domain in Y
  const double dx=lx/XMESH;	// Spatial width in X
  const double dy=ly/YMESH;	// Spatial width in Y
  const double idx=1.0/dx;
  const double idy=1.0/dy;
  const double gam=5.0/3.0;	// Specific heat ratio
  const double dr=((dx < dy)?dx:dy);
  const double cfl=0.4;		// CFL number
  const char fildir[]="dat/";	// Directory for file output
  double vfast=3.0;		// Rough estimate of vfast+vbulk for OTvortex
  double dt=cfl*dr/vfast;	// Time step

  // Variables
  double *x,*y;
  double *ro,*mx,*my,*mz,*en,*bx,*by,*bz;
}
