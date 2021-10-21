#define IFLD (6)		// Flag for initial condition. See below

void init_grid(int mpi_rank);
void init_plasma(int mpi_rank);

void init_grid(int mpi_rank)
// Define X,Y, and Z coordinates
{
  int m_xy=mpi_numx*mpi_numy;
  for (int i=0;i<nx;i++) x[i]=(i-XOFF+((mpi_rank%m_xy)%mpi_numx)*(XMESH/MNP_X))*dx;
  for (int j=0;j<ny;j++) y[j]=(j-YOFF+((mpi_rank%m_xy)/mpi_numx)*(YMESH/MNP_Y))*dy;
  for (int k=0;k<nz;k++) z[k]=(k-ZOFF+(mpi_rank/m_xy)*(ZMESH/MNP_Z))*dz;
}

void init_plasma(int mpi_rank)
// Set initial condition
{
  // Orszag-Tang vortex
  for (int k=0;k<nz;k++){
    for (int j=0;j<ny;j++){
      for (int i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;
	double vx,vy,vz,pr;

	ro[ss]=gam*gam;
	pr=gam;

#if (IFLD == 0)
	// X-Y plane
	vx=-sin(y[j]);
	vy=+sin(x[i]);
	vz=0.0;
	bx[ss]=-sin(y[j]);
	by[ss]=+sin(2*x[i]);
	bz[ss]=0.0;
#elif (IFLD ==1)
	// Y-Z plane
	vy=-sin(z[k]);
	vz=+sin(y[j]);
	vx=0.0;
	by[ss]=-sin(z[k]);
	bz[ss]=+sin(2*y[j]);
	bx[ss]=0.0;
#elif (IFLD ==2)
	// Z-X plane
	vz=-sin(x[i]);
	vx=+sin(z[k]);
	vy=0.0;
	bz[ss]=-sin(x[i]);
	bx[ss]=+sin(2*z[k]);
	by[ss]=0.0;
#elif (IFLD ==3)
	// X-Y plane (transpose)
	vx=+sin(y[j]);
	vy=-sin(x[i]);
	vz=0.0;
	bx[ss]=+sin(2*y[j]);
	by[ss]=-sin(x[i]);
	bz[ss]=0.0;
#elif (IFLD ==4)
	// Y-Z plane (transpose)
	vy=+sin(z[k]);
	vz=-sin(y[j]);
	vx=0.0;
	by[ss]=+sin(2*z[k]);
	bz[ss]=-sin(y[j]);
	bx[ss]=0.0;
#elif (IFLD ==5)
	// Z-X plane (transpose)
	vz=+sin(x[i]);
	vx=-sin(z[k]);
	vy=0.0;
	bz[ss]=+sin(2*x[i]);
	bx[ss]=-sin(z[k]);
	by[ss]=0.0;
#else
	vx=0.5*(+sin(2*z[k])-sin(y[j]));
	vy=0.5*(+sin(3*x[i])-sin(z[k]));
	vz=0.5*(+sin(4*y[j])-sin(x[i]));
	bx[ss]=0.5*(+sin(3*z[k])-sin(y[j]));
	by[ss]=0.5*(+sin(4*x[i])-sin(z[k]));
	bz[ss]=0.5*(+sin(2*y[j])-sin(x[i]));
#endif

	mx[ss]=ro[ss]*vx;
	my[ss]=ro[ss]*vy;
	mz[ss]=ro[ss]*vz;
	en[ss]=pr/(gam-1)+0.5*(ro[ss]*(vx*vx+vy*vy+vz*vz)+(bx[ss]*bx[ss]+by[ss]*by[ss]+bz[ss]*bz[ss]));
      }
    }
  }
}
