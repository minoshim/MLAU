#include "mpi.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "myfunc.h"

#define XMESH (128)		// Number of cells in X domain
#define YMESH (128)		// Number of cells in Y domain

#define MNP_X (4)		// Number of MPI processes in X
#define MNP_Y (4)		// Number of MPI processes in Y
#define MNP (MNP_X*MNP_Y)	// Total number of MPI processes

#define CFLCHECK (1)		// Flag to modify dt at every step
#define RANDOM (0)		// Flag for random noise to Vy

int dnxs[8]={0,0,0,0,0,0,0,0};
int dnys[8]={+1,+1,-1,+1,+1,+1,-1,+1};
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
  int nrec=2048;		// Number of steps for output
  int nmax=nrec*50;		// Number of maximum iteration
  const double dtrec=5.0; 	// Time step for output  
  double trec=dtrec;

  const double lx=20.0;         // Spatial domain in X
  const double ly=lx*(double)YMESH/XMESH;// Spatial domain in Y
  const double dx=lx/XMESH;	// Spatial width in X
  const double dy=ly/YMESH;	// Spatial width in Y
  const double idx=1.0/dx;
  const double idy=1.0/dy;
  const double gam=2.0;	// Specific heat ratio
  const double dr=((dx < dy)?dx:dy);
  const double cfl=0.4;		// CFL number
  const char fildir[]="dat/";	// Directory for file output

  const double beta=1e3;	// Ambient plasma beta
  const double angle_u=71.5651;	// B field angle in upper domain. 90deg: B=Bz, 0deg: B=Bx
  const double angle_l=angle_u;	// B field angle in lower domain.
  const double vamp=0.5;	// Shear velocity amplitude
  const int nmode=1;		// Number of mode for perturbation
  const double wlen=lx/nmode;
  const double lambda=1.0;	// Shear layer width
  const double s0=0.5*ly;	// Shear position
  const double ro_u=1.0;	// Density in upper domain
  const double ro_l=1.0;	// Density in lower domain
  const double b0=1.0;		// B field strength
  const double dv=0.01;		// Perturbation amplitude

  double vfast=fmode(fmin(ro_u,ro_l),beta*b0*b0,b0*b0,gam);
  double dt=cfl*dr/(vfast+vamp); // Time step

  // Variables
  double *x,*y;
  double *ro,*mx,*my,*mz,*en,*bx,*by,*bz;
}
