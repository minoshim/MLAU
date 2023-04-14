void cflcheck(double *dt, int comment)
{
  // Check CFL condition and modify time step
  int i,j,k,ss;
  double vmax=vfast;

  for (k=ZOFF;k<nz-ZOFF;k++){
    double iro,vx,vy,vz,bxc,byc,bzc,v2,b2,pr,cf,vtmp;
    for (j=YOFF;j<ny-YOFF;j++){
      for (i=XOFF;i<nx-XOFF;i++){
	ss=nx*(ny*k+j)+i;
	iro=1.0/ro[ss];
	vx=mx[ss]*iro;
	vy=my[ss]*iro;
	vz=mz[ss]*iro;
	bxc=0.5*(bx[ss]+bx[nx*(ny*k+j)+(i+1)]);
	byc=0.5*(by[ss]+by[nx*(ny*k+(j+1))+i]);
	bzc=0.5*(bz[ss]+bz[nx*(ny*(k+1)+j)+i]);
	v2=vx*vx+vy*vy+vz*vz;
	b2=bxc*bxc+byc*byc+bzc*bzc;
	pr=(gam-1)*(en[ss]-0.5*(ro[ss]*v2+b2));
	cf=fmode(ro[ss],2.0*pr,b2,gam);
	vtmp=sqrt(v2)+cf;
	if (vtmp > vmax) vmax=vtmp;
      }
    }
  }

  // MPI Allreduce
  double vmax_a;
  MPI_Allreduce(&vmax,&vmax_a,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  vmax=vmax_a;

  if (finite(vmax)) (*dt)=cfl*dr/vmax;
  if (comment){
    printf("vmax = %f, dt = %f, dx = %f\n",vmax,(*dt),dr);
    printf("CFL number = %f\n",vmax*(*dt)/dr);
  }
}
