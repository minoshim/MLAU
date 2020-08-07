#include "mpi.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "myfunc.h"

#define XMESH (2048)		// Number of cells in X domain
#define YMESH (256)		// Number of cells in Y domain
#define XOFF (4)		// Number of ghost cells in each X side.
#define YOFF (4)		// Number of ghost cells in each Y side.

#define MNP_X (8)		// Number of MPI processes in X
#define MNP_Y (2)		// Number of MPI processes in Y
#define MNP (MNP_X*MNP_Y)	// Total number of MPI processes

#define CFLCHECK (1)		// Flag to modify dt at every step
#define RANDOM (0)		// Flag for random noise to Vy
#define DIFF (1)		// Flag for diffusion
#define HALL (1)		// Flag for Hall-MHD

namespace global
{
  // Universal parameters (fixed)
  const int nx=XMESH/MNP_X+2*XOFF,ny=YMESH/MNP_Y+2*YOFF;
  const int mpi_numx=MNP_X;
  const int mpi_numy=MNP_Y;
  const double pi=4.0*atan(1.0);
  const double dtor=pi/180.;

  // Parameters
  int nrec=100;		// Number of steps for output
  int nmax=nrec*100;		// Number of maximum iteration
  const double dtrec=1.0; 	// Time step for output
  double trec=dtrec;

  const double lx=256.0;         // Spatial domain in X
  const double ly=lx*(double)YMESH/XMESH;// Spatial domain in Y
  const double xmin=-0.5*lx;		 // Minimum of x
  const double ymin=-0.5*ly;		 // Minimum of y
  const double dx=lx/XMESH;	// Spatial width in X
  const double dy=ly/YMESH;	// Spatial width in Y
  const double idx=1.0/dx;
  const double idy=1.0/dy;
  const double gam=5./3.;	// Specific heat ratio
  const double dr=((dx < dy)?dx:dy);
  const double cfl=0.2;		// CFL number
  const char fildir[]="dat/";	// Directory for file output

  // Current sheet parameters
  const double lambda=1.0;	// Current sheet thickness
  const double beta=0.2;	// Plasma beta @ lobe
  const double ro0=1.0;		// Density @ CS
  const double ro1=1.0;		// Density @ lobe
  const double b0=1.0;		// Mag field @ lobe
  const double b1=0.05;		// Mag field perturbation by Zenitani
  const double dv=0.01;		// Random noize perturbation to Vy (avaiable when RANDOM=1)

  // Diffusion coefficient
  const double eta0=DIFF*1.0*lambda/5e2; // Uniform resistivity
  const double prm=1e0;			 // Magnetic Prandtl number
  const double nu0=prm*eta0;		 // Uniform kinetic viscosity
  const double kk0=DIFF*0e-1;		 // Thermal conductivity
  const double prn=nu0/kk0;		 // Prandtl number

  // Hall parameters
  const double iner_p=0.5*lambda; // Ion inertia length
  const double eta_h=HALL*iner_p*sqrt(ro0); // Hall resistivity
  
  double vfast=fmode(ro0,beta*b0*b0,b0*b0,gam);
  double vw=eta_h*pi/dr;	// Whistler phase velocity
  double dt=cfl*dr/vfast; // Time step
  double dtw=cfl*dr/(vfast+vw);	// Time step for Hall term
  
  // Variables
  double *x,*y;
  double *ro,*mx,*my,*mz,*en,*bx,*by,*bz;
}
