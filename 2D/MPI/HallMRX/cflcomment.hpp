void cflcomment(double dx, double dt, int hflg, int mpi_rank)
// hflg: Flag to check Hall CFL
{
  // Comment CFL number
  int i,j,ss;
  double vmax=0;
  double factor=3.0;		// 3 = (Hall vel)+(Hall drift)+(whistler)
  
  for (j=YOFF;j<ny-YOFF;j++){
    double iro,vx,vy,vz,bxc,byc,bzc,v2,b2,pr,cf,vtmp;
    for (i=XOFF;i<nx-XOFF;i++){
      ss=nx*j+i;
      iro=1.0/ro[ss];
      vx=mx[ss]*iro;
      vy=my[ss]*iro;
      vz=mz[ss]*iro;
      bxc=0.5*(bx[ss]+bx[nx*j+(i+1)]);
      byc=0.5*(by[ss]+by[nx*(j+1)+i]);
      bzc=bz[ss];
      v2=vx*vx+vy*vy+vz*vz;
      b2=bxc*bxc+byc*byc+bzc*bzc;
      pr=(gam-1)*(en[ss]-0.5*(ro[ss]*v2+b2));
      cf=fmode(ro[ss],2.0*pr,b2,gam);
      vtmp=sqrt(v2)+cf+hflg*vw*sqrt(b2)*iro*factor;
      if (vtmp > vmax) vmax=vtmp;
    }
  }

  // MPI Allreduce
  double vmax_a;
  MPI_Allreduce(&vmax,&vmax_a,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  vmax=vmax_a;

  if (mpi_rank ==0){
    printf("vmax = %f, dt = %f, dx = %f\n",vmax,dt,dx);
    printf("CFL number = %f\n",vmax*dt/dx);
    if (hflg == 0) printf("CFL condition for diffusion is %f\n ",fmax(eta0,nu0)*dt/(dx*dx));
  }
}
