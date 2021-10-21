void init_grid(int mpi_rank);
void init_plasma(int mpi_rank);

void init_grid(int mpi_rank)
// Define X and Y coordinates
{
  for (int i=0;i<nx;i++) x[i]=(i-XOFF+0.5+(mpi_rank%mpi_numx)*(XMESH/MNP_X))*dx;
  for (int j=0;j<ny;j++) y[j]=(j-YOFF+0.5+(mpi_rank/mpi_numx)*(YMESH/MNP_Y))*dy;
}

void init_plasma(int mpi_rank)
// Set initial condition
{
  // KHI
  int i,j,ss;
  double dvy[nx];
  double stim;			// Time at node 0, used for seed of rand_noise
  unsigned seed;

  if (mpi_rank == 0) stim=MPI_Wtime();
  MPI_Bcast(&stim,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  seed=(unsigned)(stim*(1+(mpi_rank%mpi_numx)));

  for (i=0;i<nx;i++){
    dvy[i]=0.0;
#if (RANDOM)
    double dvpara[2]={0,dv};
    dvy[i]+=rand_noise(dvpara,seed); // Multiple mode perturbation
#else
    dvy[i]+=dv*sin(2*pi*x[i]/wlen); // Single mode perturbation
#endif    
  }

  for (j=0;j<ny;j++){
    double angle=0.5*((angle_u+angle_l)+(angle_u-angle_l)*tanh((y[j]-s0)/lambda));
    for (i=0;i<nx;i++){
      double vx,vy,vz,pr,bxc,byc,bzc,b2;
      double xm,ym;

      ss=nx*j+i;
      xm=x[i]-0.5*dx;
      ym=y[j]-0.5*dy;

      ro[ss]=0.5*((ro_u+ro_l)+(ro_u-ro_l)*tanh((y[j]-s0)/lambda));

      vx=+vamp*tanh((y[j]-s0)/lambda);
      vy=dvy[i]*exp(-((y[j]-s0)*(y[j]-s0))/(4*lambda*lambda));
      vz=0.0;

      if (angle == 90){
	bx[ss]=bxc=0.0;
	by[ss]=byc=0.0;
	bz[ss]=bzc=b0;
      } else{
	bx[ss]=bxc=b0*cos(angle*dtor);
	by[ss]=byc=0.0;
	bz[ss]=bzc=b0*sin(angle*dtor);
      }
      b2=bxc*bxc+byc*byc+bzc*bzc;
      
      pr=0.5*beta*b2;

      mx[ss]=ro[ss]*vx;
      my[ss]=ro[ss]*vy;
      mz[ss]=ro[ss]*vz;
      en[ss]=pr/(gam-1)+0.5*(ro[ss]*(vx*vx+vy*vy+vz*vz)+b2);
    }
  }
}
