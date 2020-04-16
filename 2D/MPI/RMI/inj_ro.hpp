void inj_ro(double *ro, double *mx, double *my, double *mz,
	    double *en, double *bx, double *by, double *bz,
	    double ro_3, double dro3, double vy_1, double bx_1, double by_1, double bz_1, double pr_1,
	    int nx, int ny, int xoff, int yoff,
	    int mpi_rank, int mpi_numx, int mpi_numy)
// Set density @ upper boundary for injection
// ro_3: CD density @ upper boundary
// dro3: CD density perturbation @ upper boundary
// vy_1, bx_1, by_1, bz_1, pr_1: upstream variables
{
   
  int i,j,ss;
  static int irflg=0;
  
  if (irflg == 0){
    double stim;			// Time at node 0, used for seed of random
    unsigned seed;
    if (mpi_rank == 0) stim=MPI_Wtime();
    MPI_Bcast(&stim,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    seed=(unsigned)(stim*(1+mpi_rank));
    srandom(seed);
    irflg=1;
  }

  if ((mpi_rank / mpi_numx) == (mpi_numy-1)){
    for (j=ny-YOFF;j<ny;j++){
      int ctyflg=(j == (ny-YOFF));
      for (i=0;i<nx;i++){
	double vx,vy,vz,pr,b2;
	double bxc,byc,bzc;	// Centered quantities
	ss=nx*j+i;
	
	ro[ss]=ro_3;
#if (RANDOM)
	ro[ss]+=dro3*((double)random()/RAND_MAX-0.5)*2.0;
#endif
	vx=0;
	vy=vy_1;
	vy-=vref;
	vz=0;

	bx[ss]=bxc=bx_1;
	by[ss]=ctyflg?by[ss]:by_1; // Due to CT spacing
	byc=by_1;
	bz[ss]=bzc=bz_1;
	b2=bxc*bxc+byc*byc+bzc*bzc;

	pr=pr_1;

	mx[ss]=ro[ss]*vx;
	my[ss]=ro[ss]*vy;
	mz[ss]=ro[ss]*vz;
	en[ss]=pr/(gam-1)+0.5*(ro[ss]*(vx*vx+vy*vy+vz*vz)+b2);
      }
    }
  }
  
}
