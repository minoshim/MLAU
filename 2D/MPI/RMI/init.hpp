void init_grid(int mpi_rank);
void init_plasma(int mpi_rank);

void init_grid(int mpi_rank)
// Define X and Y coordinates
{
  for (int i=0;i<nx;i++) x[i]=(i-xoff+0.5+(mpi_rank%mpi_numx)*(XMESH/MNP_X))*dx+xmin;
  for (int j=0;j<ny;j++) y[j]=(j-yoff+0.5+(mpi_rank/mpi_numx)*(YMESH/MNP_Y))*dy+ymin;
}

void init_plasma(int mpi_rank)
// Set initial condition
{
  int i,j,ss;
  double dw=0.5*dy;
  
  // R-H relation
  double a=2.0*(2.0-gam)/beta;
  double b=gam*((gam-1.0)*ma_u*ma_u+2.0/beta+2.0);
  double c=-gam*(gam+1.0)*ma_u*ma_u;
  // Compression ratio
  double r=(-b+sqrt(b*b-4.0*a*c))/(2.0*a);
  // Pr2/Pr1
  double rr=gam*ma_u*ma_u*(1.0-1.0/r)-(r*r-1.0)/beta+1.0;
  // Donwstream paramters
  double ro_2=ro_1*r;
  double vy_2=vy_1/r;
  double bx_2=bx_1*r;
  double bz_2=bz_1*r;
  double pr_2=pr_1*rr;

  // Messaage of initial condition
  if (mpi_rank ==0){
    puts("Initial condition in shock rest frame");
    printf("U: %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",ro_1,0.,vy_1,0.,bx_1,by_1,bz_1,pr_1);
    printf("D: %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",ro_2,0.,vy_2,0.,bx_2,by_1,bz_2,pr_2);
    printf("Mach = %.9f, beta = %.9f\n",ma_u,beta);
  }

  // Random seed for contact discon.
  double stim;			// Time at node 0, used for seed of rand_noise
  unsigned seed;
  if (mpi_rank == 0) stim=MPI_Wtime();
  MPI_Bcast(&stim,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  seed=(unsigned)(stim*(1+mpi_rank));

  // Position of contact discon.
  double ypos[nx];
  for (i=0;i<nx;i++){
    ypos[i]=1.0+psi*cos(2*pi*x[i]/lambda); // Single mode
  }
#if (RANDOM)
  // Multi-mode
  double para[2]={pi,pi};
  for (int m=2;m<=mmax;m++){
    double xphase;
    if (mpi_rank == 0) xphase=rand_noise(para,seed);
    MPI_Bcast(&xphase,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (i=0;i<nx;i++){
      // ypos[i]+=(psi/m)*cos(2*pi*m*(x[i]-xphase)/lambda);
      ypos[i]+=(psi/(m*m))*cos(2*pi*m*(x[i]-xphase)/lambda);
    }
  }
#endif

  for (j=0;j<ny;j++){
    double sfunc;    
    sfunc=0.5*(1.0+tanh(y[j]/dw)); // 0 (y<0), 1 (y>0)
    for (i=0;i<nx;i++){
      double vx,vy,vz,pr,b2;
      double bxc,byc,bzc;	// Centered quantities
      ss=nx*j+i;
      ro[ss]=ro_2+(ro_1-ro_2)*sfunc;

      vx=0;      
      vy=vy_2+(vy_1-vy_2)*sfunc;
      vy-=vref;
      vz=0;

      bx[ss]=bxc=bx_2+(bx_1-bx_2)*sfunc;
      by[ss]=byc=by_1;
      bz[ss]=bzc=bz_2+(bz_1-bz_2)*sfunc;
      b2=bxc*bxc+byc*byc+bzc*bzc;
      
      pr=pr_2+(pr_1-pr_2)*sfunc;
      
      // Contact discon
      double rocd=ro_3;
#if (RANDOM)
      double rpara[2]={ro_3,dro3};
      rocd=rand_noise(rpara,seed);
#endif
      ro[ss]+=(rocd-ro_1)*0.5*(1.0+tanh((y[j]-ypos[i])/dw));

      mx[ss]=ro[ss]*vx;
      my[ss]=ro[ss]*vy;
      mz[ss]=ro[ss]*vz;
      en[ss]=pr/(gam-1)+0.5*(ro[ss]*(vx*vx+vy*vy+vz*vz)+b2);
    }
  }

}
