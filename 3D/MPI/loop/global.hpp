#include "mpi.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "myfunc.h"

#define XMESH (100)		// Number of cells in X domain
#define YMESH (100)		// Number of cells in Y domain
#define ZMESH (200)		// Number of cells in Z domain

#define MNP_X (2)		// Number of MPI processes in X
#define MNP_Y (2)		// Number of MPI processes in Y
#define MNP_Z (2)		// Number of MPI processes in Z
#define MNP (MNP_X*MNP_Y*MNP_Z)	// Total number of MPI processes

#define CFLCHECK (0)		// Flag to modify dt at every step

int dnxs[8]={0,0,0,0,0,0,0,0};
int dnys[8]={0,0,0,0,0,0,0,0};
int dnzs[8]={0,0,0,0,0,0,0,0};
// Boundary condition flag for ro,mx,my,mz,en,bx,by,bz (be sure of variable order)
// 0=Periodic, +1=Neumann, -1=Dirichlet
int stxs[8]={0,0,0,0,0,1,0,0};
int stys[8]={0,0,0,0,0,0,1,0};
int stzs[8]={0,0,0,0,0,0,0,1};
// Staggered grid flag for ro,mx,my,mz,en,bx,by,bz (be sure of variable order)
// Do NOT change

namespace global
{
  // Universal parameters (fixed)
  const int xoff=4,yoff=xoff,zoff=xoff;	// Number of ghost cells in each side
  const int nx=XMESH/MNP_X+2*xoff,ny=YMESH/MNP_Y+2*yoff,nz=ZMESH/MNP_Z+2*zoff;
  const int mpi_numx=MNP_X;
  const int mpi_numy=MNP_Y;
  const int mpi_numz=MNP_Z;
  const double pi=4.0*atan(1.0);
  const double dtor=pi/180.;

  // Parameters
  int nrec=200;			// Number of steps for output
  int nmax=nrec*5;		// Number of maximum iteration
  const double dtrec=0.2; 	// Time step for output  
  double trec=dtrec;

  const double lx=1.0;	// Spatial domain in X
  const double ly=1.0;	// Spatial domain in Y
  const double lz=2.0;	// Spatial domain in Z
  const double xmin=-0.5*lx;	// Leftmost x value
  const double ymin=-0.5*ly;	// Leftmost y value
  const double zmin=-0.5*lz;	// Leftmost z value
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

  // looe advection of Lee13
  const double ro0=1.0;
  const double pr0=1e0;
  const double v_0=1.0;
  const double a0=1e-3;
  const double rad=0.3;
  const double tilt=atan(0.5);	// Loop tilt around y-axis
  const double theta=0.01; // Advection angle (radian)
  const double vx0=v_0*cos(theta);
  const double vy0=v_0*sin(theta);
  const double vz0=v_0*2.0;
  
  double vfast=fmax(fmax(v_0,vz0),sqrt(pr0/ro0));
  double dt=cfl*dr/vfast;	// Time step

  // Variables
  double *x,*y,*z;
  double *ro,*mx,*my,*mz,*en,*bx,*by,*bz;
}
