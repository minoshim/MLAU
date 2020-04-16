#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "myfunc.h"

#define XMESH (800)		// Number of cells in computational domain
#define XOFF (4)		// Number of ghost cells in each side. fix to be 4.

#define NUM (1)			// Select initial condition
// NUM = 1: Dai & Woodward 1994 (Miyoshi & Kusano 2005, Fig. 5)
// NUM = 2: Brio & Wu 1988 (Miyoshi & Kusano 2005, Fig. 8)
// NUM = 3: Slow switch-off shock (Miyoshi & Kusano 2005, Fig. 9)
// NUM = 4: Slow switch-off rarefaction (Miyoshi & Kusano 2005, Fig. 10)
// NUM = 5: Super-fast expansion (Miyoshi & Kusano 2005, Fig. 11)

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
  double bx;			// Normal magnetic field strength
  double gamma=5./3.;	// Specific heat ratio
  const double cfl=0.4;		// CFL number
  const char fildir[]="dat/";	// Directory for file output
  char filname[100];
  FILE *outfil;
  double x[nx],ro[nx],mx[nx],my[nx],mz[nx],en[nx],by[nx],bz[nx];

  // Set normal magnetic field
#if (NUM == 2)
  // BrioWu
  bx=0.75;
  gamma=2.0;
#elif (NUM == 3)
  // Slow switch-off shock
  bx=1.0;
#elif (NUM == 4)
  // Slow switch-off rarefaction
  bx=1.0;
#elif (NUM == 5)
  // Super-fast expansion
  bx=0.0;
#else  // Default
  // DaiWoodward
  bx=2.0/sqrt(4*pi);
#endif

  /* Initialize */
  for (i=0;i<nx;i++){
    double vx,vy,vz,pr;
    x[i]=(i-XOFF+0.5)*dx-0.5;
    char flag=(x[i] <= 0);
    
#if (NUM == 2)
    // BrioWu
    ro[i]=(flag)?(1.00):(0.125);
    vx=0;
    vy=0;
    vz=0;
    by[i]=(flag)?(1):(-1);
    bz[i]=0;
    pr=(flag)?(1.0):(0.1);
#elif (NUM ==3)
  // Slow switch-off shock
    ro[i]=(flag)?(1.368):(1.0);
    vx=(flag)?(0.269):(0.0);
    vy=(flag)?(1.0):(0.0);
    vz=(flag)?(0.0):(0.0);
    by[i]=(flag)?(0.0):(1.0);
    bz[i]=0.0;
    pr=(flag)?(1.769):(1.0);
#elif (NUM == 4)
  // Slow switch-off rarefaction
    ro[i]=(flag)?(1.0):(0.2);
    vx=(flag)?(0.0):(1.186);
    vy=(flag)?(0.0):(2.967);
    vz=(flag)?(0.0):(0.0);
    by[i]=(flag)?(0.0):(1.6405);
    bz[i]=0.0;
    pr=(flag)?(2.0):(0.1368);
#elif (NUM == 5)    
  // Super-fast expansion
    ro[i]=1.0;
    vx=(flag)?(-3.1):(3.1);
    vy=0.0;
    vz=0.0;
    by[i]=0.5;
    bz[i]=0.0;
    pr=0.45;
#else  // Default
    // DaiWoodward
    ro[i]=(flag)?(1.08):(1.0);
    vx=(flag)?(1.20):(0.0);
    vy=(flag)?(0.01):(0.0);
    vz=(flag)?(0.50):(0.0);
    by[i]=(flag)?(3.6/sqrt(4*pi)):(4.0/sqrt(4*pi));
    bz[i]=2.0/sqrt(4*pi);
    pr=(flag)?(0.95):(1.0);
#endif

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
  int nmax=nrec*50;		// Number of maximum iteration
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
