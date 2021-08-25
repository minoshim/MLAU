#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "myfunc.h"

#define XMESH (256)		// Number of cells in computational domain
#define XOFF (4)		// Number of ghost cells in each side. fix to be 4.

int main(void)
{
  int i,n=0,cnt=0;
  const int nx=XMESH+2*XOFF;	// Number of cells in whole domain (including offset)
  const double pi=4.0*atan(1.0);  
  double tim=0.0;		// Time
  double dtrec=0.01;		// Time step for output
  double dt;			// Time step for calculation
  const double lx=1.0;		// Spatial domain size
  const double dx=lx/XMESH;	// Spatial width
  const double bx=1.0;		// Normal magnetic field strength
  const double gamma=5./3.;	// Specific heat ratio
  const double beta=10.0;	// Plasma beta
  const double ro0=1.0;		// Ambient density
  const double pr0=0.5*beta*bx*bx; // Ambient pressure
  const double cf0=sqrt((bx*bx+gamma*pr0)/ro0); // Magnetosonic speed
  const double vx0=0.0*cf0;	// Ambient normal velocity
  const double cfl=0.4;		// CFL number
  const char fildir[]="dat/";	// Directory for file output
  char filname[100];
  FILE *outfil;
  double x[nx],ro[nx],mx[nx],my[nx],mz[nx],en[nx],by[nx],bz[nx];

  printf("vx0 = %f cf0 = %f\n",vx0,cf0);
  
  /* Initialize */
  unsigned seed=10;
  // seed=(unsigned)time(NULL);
  for (i=0;i<nx;i++){
    double vx,vy,vz,pr;
    x[i]=(i-XOFF+0.5)*dx;

    double para1[]={0,0.01},para2[]={1,0.01};
    ro[i]=ro0;
    vx=vx0;
    vy=0.0;
    vz=0.0;
    by[i]=rand_noise(para1,seed);
    bz[i]=rand_noise(para1,seed);
    pr=pr0*rand_noise(para2,seed);
    
    mx[i]=ro[i]*vx;
    my[i]=ro[i]*vy;
    mz[i]=ro[i]*vz;
    en[i]=pr/(gamma-1)+0.5*(ro[i]*(vx*vx+vy*vy+vz*vz)+(bx*bx+by[i]*by[i]+bz[i]*bz[i]));
  }

  // Calculate time step
  double vmax=1.0;
  for (i=0;i<nx;i++){
    double iro,vx,vy,vz,v2,b2,pr,cf,vtmp;
    iro=1.0/ro[i];
    vx=mx[i]*iro;
    vy=my[i]*iro;
    vz=mz[i]*iro;
    v2=vx*vx+vy*vy+vz*vz;
    b2=bx*bx+by[i]*by[i]+bz[i]*bz[i];
    pr=(gamma-1)*(en[i]-0.5*(ro[i]*v2+b2));
    cf=sqrt((gamma*pr+b2)*iro);
    vtmp=fabs(vx)+cf;
    if (vtmp > vmax) vmax=vtmp;
  }
  dt=cfl*dx/vmax;
  int nrec=(int)(dtrec/dt+0.5);	// Number of iterations for output
  int nmax=nrec*4095;		// Number of maximum iteration
  printf("Data output every %d steps (%f duration) \n",nrec,dtrec);

  /* Output */
  sprintf(filname,"%s/x.dat",fildir);
  outfil=fopen(filname,"w");
  for (i=0;i<nx;i++) fprintf(outfil,"%.12f\n",x[i]);
  fclose(outfil);
  sprintf(filname,"%s/params.dat",fildir);  
  outfil=fopen(filname,"w");
  fprintf(outfil,"%.12f %.12f\n",bx,gamma);
  fclose(outfil);
  sprintf(filname,"%s/xoff.dat",fildir);  
  outfil=fopen(filname,"w");
  fprintf(outfil,"%d\n",XOFF);
  fclose(outfil);
  sprintf(filname,"%s/t.dat",fildir);  
  outfil=fopen(filname,"w");
  fprintf(outfil,"%.12f\n",tim);
  fclose(outfil);
  sprintf(filname,"%s/outdat_%05d.dat",fildir,cnt);  
  outfil=fopen(filname,"wb");
  fwrite(&ro[0],sizeof(double),nx,outfil);
  fwrite(&mx[0],sizeof(double),nx,outfil);
  fwrite(&my[0],sizeof(double),nx,outfil);
  fwrite(&mz[0],sizeof(double),nx,outfil);
  fwrite(&en[0],sizeof(double),nx,outfil);
  fwrite(&by[0],sizeof(double),nx,outfil);
  fwrite(&bz[0],sizeof(double),nx,outfil);
  fclose(outfil);

  /* Time integration */
  while(n++ < nmax){
    tim+=dt;

    // Boundary condition
    prdc_1d(ro,nx,XOFF);
    prdc_1d(mx,nx,XOFF);
    prdc_1d(my,nx,XOFF);
    prdc_1d(mz,nx,XOFF);
    prdc_1d(en,nx,XOFF);
    prdc_1d(by,nx,XOFF);
    prdc_1d(bz,nx,XOFF);

    // MHD update
    mhd_fd4c_1d(&ro[0],&mx[0],&my[0],&mz[0],&en[0],&by[0],&bz[0],
		bx,dt,dx,nx,XOFF,gamma);

    /* Output */
    if ((n % nrec) == 0){
      cnt++;
      sprintf(filname,"%s/t.dat",fildir);  
      outfil=fopen(filname,"a");
      fprintf(outfil,"%.12f\n",tim);
      fclose(outfil);
      sprintf(filname,"%s/outdat_%05d.dat",fildir,cnt);  
      outfil=fopen(filname,"wb");
      fwrite(&ro[0],sizeof(double),nx,outfil);
      fwrite(&mx[0],sizeof(double),nx,outfil);
      fwrite(&my[0],sizeof(double),nx,outfil);
      fwrite(&mz[0],sizeof(double),nx,outfil);
      fwrite(&en[0],sizeof(double),nx,outfil);
      fwrite(&by[0],sizeof(double),nx,outfil);
      fwrite(&bz[0],sizeof(double),nx,outfil);
      fclose(outfil);
    }
  }
  return 0;
}
