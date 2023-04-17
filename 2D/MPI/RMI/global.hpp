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
#define YMESH (XMESH*80)		// Number of cells in Y domain

#define MNP_X (2)		// Number of MPI processes in X
#define MNP_Y (8)		// Number of MPI processes in Y
#define MNP (MNP_X*MNP_Y)	// Total number of MPI processes

#define CFLCHECK (1)		// Flag to modify dt at every step
#define RANDOM (0)		// Flag for random noise

int dnxs[8]={0,0,0,0,0,0,0,0};
int dnys[8]={+1,+1,+1,+1,+1,+1,+1,+1};
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
  int nrec=400;		// Number of steps for output
  int nmax=nrec*50;		// Number of maximum iteration
  const double dtrec=0.5; 	// Time step for output  
  double trec=dtrec;

  const double lx=1.0;         // Spatial domain in X
  const double ly=lx*(double)YMESH/XMESH;// Spatial domain in Y
  const double xmin=0.0;		 // Minimum of X
  const double ymin=-0.5*ly;		 // Minimum of Y
  const double dx=lx/XMESH;	// Spatial width in X
  const double dy=ly/YMESH;	// Spatial width in Y
  const double idx=1.0/dx;
  const double idy=1.0/dy;
  const double gam=5./3.;	// Specific heat ratio
  const double dr=((dx < dy)?dx:dy);
  const double cfl=0.4;		// CFL number
  const char fildir[]="dat/";	// Directory for file output

  // Paramters for RMI
  //// Shock upstream paramters
  const double ro_1=+1.0;	// Density
  const double vy_1=-1.0;	// Velocity
  const double ma_u=100.0;	// Sound Mach number
  const double beta=1e5;	// Beta
  const double pr_1=(ro_1*vy_1*vy_1)/(gam*ma_u*ma_u); // Pressure
  const double b0_1=sqrt(2.0*pr_1/beta);	      // Magnetic field
  const double bx_1=b0_1;
  const double by_1=0.0;
  const double bz_1=0.0;
  const double vref=-0.625;	// Velocity of reference. 0 = shock rest frame
  //// Contact Discon parameters
  const double ro_3=ro_1*1e1;	// Density jump
  const double lambda=lx;	// Wavelength
  const double psi=0.1*lambda;	// Amplitude
  const double mmax=4;		// Number of modes (Available when RAMDOM=1)
  const double dro3=ro_3*0.1;	// Density perturbation (Available when RAMDOM=1)

  const double vfast=fmode(ro_1,2.0*pr_1,b0_1*b0_1,gam);
  double dt=cfl*dr/vfast; // Time step

  // Variables
  double *x,*y;
  double *ro,*mx,*my,*mz,*en,*bx,*by,*bz;
}
