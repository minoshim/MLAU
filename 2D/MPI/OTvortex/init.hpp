void init_grid(int mpi_rank);
void init_plasma(int mpi_rank);

void init_grid(int mpi_rank)
// Define X and Y coordinates
{
  for (int i=0;i<nx;i++) x[i]=(i-XOFF+(mpi_rank%mpi_numx)*(XMESH/MNP_X))*dx;
  for (int j=0;j<ny;j++) y[j]=(j-YOFF+(mpi_rank/mpi_numx)*(YMESH/MNP_Y))*dy;
}

void init_plasma(int mpi_rank)
// Set initial condition
{
  // Orszag-Tang vortex
  for (int j=0;j<ny;j++){
    for (int i=0;i<nx;i++){
      int ss=nx*j+i;
      double vx,vy,vz,pr;
      ro[ss]=gam*gam;
      vx=-sin(y[j]);
      vy=+sin(x[i]);
      vz=0.0;
      mx[ss]=ro[ss]*vx;
      my[ss]=ro[ss]*vy;
      mz[ss]=ro[ss]*vz;
      pr=gam;
      bx[ss]=-sin(y[j]);
      by[ss]=+sin(2*x[i]);
      bz[ss]=0.0;
      en[ss]=pr/(gam-1)+0.5*(ro[ss]*(vx*vx+vy*vy+vz*vz)+(bx[ss]*bx[ss]+by[ss]*by[ss]+bz[ss]*bz[ss]));
    }
  }
}
