void cflcheck(double dr, double *dt, int hflg)
// hflg: Flag to check Hall CFL
{
  // Check CFL condition and modify time step
  int i,j,ss;
  double vmax=vfast+hflg*vw;
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

  if (finite(vmax)) (*dt)=cfl*dr/vmax;
}
