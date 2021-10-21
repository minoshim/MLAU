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
    x[i]=(i-XOFF+0.5)*dx-0.5*lx;
    fprintf(outfil,"%.12f\n",x[i]);
  }
  fclose(outfil);

  sprintf(filname,"%s/y.dat",fildir);
  outfil=fopen(filname,"w");
  for (int j=0;j<ny;j++){
    y[j]=(j-YOFF+0.5)*dy-0.5*ly;
    fprintf(outfil,"%.12f\n",y[j]);
  }
  fclose(outfil);

}

void init_plasma()
// Set initial condition
{
  // Blast wave
  for (int j=0;j<ny;j++){
    for (int i=0;i<nx;i++){
      int ss=nx*j+i;
      double rr,vx,vy,vz,pr;
      rr=sqrt(x[i]*x[i]+y[j]*y[j]);
      
      ro[ss]=(rr <= r_0)?ro1:ro0;
      vx=0.0;
      vy=0.0;
      vz=0.0;
      pr=(rr <= r_0)?pr1:pr0;
      
      bx[ss]=b_0*cos(ban*dtor);
      by[ss]=b_0*sin(ban*dtor);
      bz[ss]=0.0;

      mx[ss]=ro[ss]*vx;
      my[ss]=ro[ss]*vy;
      mz[ss]=ro[ss]*vz;
      en[ss]=pr/(gam-1)+0.5*(ro[ss]*(vx*vx+vy*vy+vz*vz)+(bx[ss]*bx[ss]+by[ss]*by[ss]+bz[ss]*bz[ss]));
    }
  }
}
