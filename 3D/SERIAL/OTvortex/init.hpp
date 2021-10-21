#define IFLD (0)		// Flag for initial condition. See init_plasma

void init_grid();
void init_plasma();

void init_grid()
// Define X,Y and Z coordinates
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

  sprintf(filname,"%s/z.dat",fildir);
  outfil=fopen(filname,"w");
  for (int k=0;k<nz;k++){
    // z[k]=(k-ZOFF+0.5)*dz;
    z[k]=(k-ZOFF)*dz;
    fprintf(outfil,"%.12f\n",z[k]);
  }
  fclose(outfil);

}

void init_plasma()
// Set initial condition
{
  // Orszag-Tang vortex
  for (int k=0;k<nz;k++){
    for (int j=0;j<ny;j++){
      for (int i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;
	double vx,vy,vz,pr;

	ro[ss]=gam*gam;
	pr=gam;

#if (IFLD == 0)
	// X-Y plane
	vx=-sin(y[j]);
	vy=+sin(x[i]);
	vz=0.0;
	bx[ss]=-sin(y[j]);
	by[ss]=+sin(2*x[i]);
	bz[ss]=0.0;
#elif (IFLD ==1)
	// Y-Z plane
	vy=-sin(z[k]);
	vz=+sin(y[j]);
	vx=0.0;
	by[ss]=-sin(z[k]);
	bz[ss]=+sin(2*y[j]);
	bx[ss]=0.0;
#elif (IFLD ==2)
	// Z-X plane
	vz=-sin(x[i]);
	vx=+sin(z[k]);
	vy=0.0;
	bz[ss]=-sin(x[i]);
	bx[ss]=+sin(2*z[k]);
	by[ss]=0.0;
#elif (IFLD ==3)
	// X-Y plane (transpose)
	vx=+sin(y[j]);
	vy=-sin(x[i]);
	vz=0.0;
	bx[ss]=+sin(2*y[j]);
	by[ss]=-sin(x[i]);
	bz[ss]=0.0;
#elif (IFLD ==4)
	// Y-Z plane (transpose)
	vy=+sin(z[k]);
	vz=-sin(y[j]);
	vx=0.0;
	by[ss]=+sin(2*z[k]);
	bz[ss]=-sin(y[j]);
	bx[ss]=0.0;
#elif (IFLD ==5)
	// Z-X plane (transpose)
	vz=+sin(x[i]);
	vx=-sin(z[k]);
	vy=0.0;
	bz[ss]=+sin(2*x[i]);
	bx[ss]=-sin(z[k]);
	by[ss]=0.0;
#else
	vx=0.5*(+sin(2*z[k])-sin(y[j]));
	vy=0.5*(+sin(3*x[i])-sin(z[k]));
	vz=0.5*(+sin(4*y[j])-sin(x[i]));
	bx[ss]=0.5*(+sin(3*z[k])-sin(y[j]));
	by[ss]=0.5*(+sin(4*x[i])-sin(z[k]));
	bz[ss]=0.5*(+sin(2*y[j])-sin(x[i]));
#endif

	mx[ss]=ro[ss]*vx;
	my[ss]=ro[ss]*vy;
	mz[ss]=ro[ss]*vz;
	en[ss]=pr/(gam-1)+0.5*(ro[ss]*(vx*vx+vy*vy+vz*vz)+(bx[ss]*bx[ss]+by[ss]*by[ss]+bz[ss]*bz[ss]));
      }
    }
  }
}
