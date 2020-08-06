void init_plasma()
// Set initial condition
{
  // Orszag-Tang vortex
  for (int j=0;j<ny;j++){
    for (int i=0;i<nx;i++){
      int ss=nx*j+i;
      double vx,vy,vz,pr;
      ro[ss]=ro0;
      vx=-v0*sin(2.0*pi*y[j]);
      vy=+v0*sin(2.0*pi*x[i]);
      vz=0.0;
      mx[ss]=ro[ss]*vx;
      my[ss]=ro[ss]*vy;
      mz[ss]=ro[ss]*vz;
      pr=pr0;
      bx[ss]=-b0*sin(2.0*pi*y[j]);
      by[ss]=+b0*sin(4.0*pi*x[i]);
      bz[ss]=0.0;
      en[ss]=pr/(gam-1)+0.5*(ro[ss]*(vx*vx+vy*vy+vz*vz)+(bx[ss]*bx[ss]+by[ss]*by[ss]+bz[ss]*bz[ss]));
    }
  }
}
