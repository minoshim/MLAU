void init_grid();
void init_plasma();

void init_grid()
// Define X and Y coordinates
{
  char filname[100];
  FILE *outfil;

  sprintf(filname,"%s/x.dat",fildir);
  outfil=fopen(filname,"w");
  for (int i=0;i<nx;i++){
    // x[i]=(i-XOFF+0.5)*dx;
    x[i]=(i-XOFF)*dx;
    fprintf(outfil,"%.12f\n",x[i]);
  }
  fclose(outfil);

  sprintf(filname,"%s/y.dat",fildir);
  outfil=fopen(filname,"w");
  for (int j=0;j<ny;j++){
    // y[j]=(j-YOFF+0.5)*dy;
    y[j]=(j-YOFF)*dy;
    fprintf(outfil,"%.12f\n",y[j]);
  }
  fclose(outfil);

}

void init_plasma()
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
