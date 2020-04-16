/* Copyright 2020 Takashi Minoshima */

/* This file is part of MLAU. */

/* MLAU is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* MLAU is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with MLAU.  If not, see <https://www.gnu.org/licenses/>. */

#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "interp.h"
#include "mhd_flux.h"
#include "mhd_eigen.h"

#define ODR (1)			/* Spatial order (1,2,3,4) */
#define R_K (2)			/* Runge-Kutta order (1,2,3). Never set >3 */

void openbc1d_mhd(double *ro, double *mx, double *my, double *mz,
		  double *en, double *by, double *bz,
		  double bx,
		  int nx, int xoff, double gamma,
		  int istt, int iend);

inline void calc_flr(const double *f, double *fl, double *fr)
/* Interpolation to cell face */
{
#if (ODR == 1)
  *fl=f[2]; *fr=f[2];
#elif (ODR == 2)
  muscl_mm_cal_flr(f[1],f[2],f[3],fl,fr); /* MinMod */
  /* muscl_mc_cal_flr(f[1],f[2],f[3],fl,fr); /\* MC *\/ */
#elif (ODR == 3)
  wcns3_cal_flr(f[1],f[2],f[3],fl,fr);
#else
  wcns4_cal_flr(f[0],f[1],f[2],f[3],f[4],fl,fr);
#endif
}

inline double calc_df(const double *f)
/* Evaluate 1st difference */
{
#if (ODR <= 2)
  return(f[2]-f[1]);
#else
  return((27.0*(f[2]-f[1])-(f[3]-f[0]))/24.0);
#endif
}

void mhd_fd4c_1d(double *ro, double *mx, double *my, double *mz,
		 double *en, double *by, double *bz,
		 double bx, double dt, double dx,
		 int nx, int xoff, double gamma)
/* 1D MHD simulation  */
/* ro: density */
/* m?: moment */
/* en: energy normalized by B^2/(4 pi) */
/* b?: magnetic field */
/* gamma: specific heat ratio */
{
  int i,ss,rk;
  const int ns=5;		/* number of grids in the stencil */
  const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const double dtdx=dt/dx;
  const double bx2=bx*bx;
  double *ut,*ul,*ur,*fx;
  double *vx,*vy,*vz,*pr; 
  void (*func_lr)(const double*, double*, double*)=&calc_flr;
  double (*func_df)(const double*)=&calc_df;
  void (*func_flux)(double, double, double, double, double, double, double, 
		    double, double, double, double, double, double, double, 
		    double, double, double*,
		    double*, double*, double*, double*, double*, double*, double*);
  func_flux=&calc_flux_mlau;
  /* func_flux=&calc_flux_hlld; */
  /* func_flux=&calc_flux_roe; */

  ut=(double*)malloc(sizeof(double)*7*nx);
  ul=(double*)malloc(sizeof(double)*7*nx);
  ur=(double*)malloc(sizeof(double)*7*nx);
  fx=(double*)malloc(sizeof(double)*7*nx);
  vx=(double*)malloc(sizeof(double)*nx);
  vy=(double*)malloc(sizeof(double)*nx);
  vz=(double*)malloc(sizeof(double)*nx);
  pr=(double*)malloc(sizeof(double)*nx);

  /* Copy */
  for (i=0;i<nx;i++){
    ss=i;
    ut[7*ss+0]=ro[ss];
    ut[7*ss+1]=mx[ss];
    ut[7*ss+2]=my[ss];
    ut[7*ss+3]=mz[ss];
    ut[7*ss+4]=by[ss];
    ut[7*ss+5]=bz[ss];
    ut[7*ss+6]=en[ss];
  }

  for (rk=0;rk<R_K;rk++){
#ifdef _OPENMP
#pragma omp parallel private(i,ss)
#endif
    {

      /* Primitive variable at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
#pragma simd
      for (i=0;i<nx;i++){
	double iro,v2,b2;
	ss=i;
	iro=1.0/ro[ss];
	vx[ss]=mx[ss]*iro;
	vy[ss]=my[ss]*iro;
	vz[ss]=mz[ss]*iro;
	v2=vx[ss]*vx[ss]+vy[ss]*vy[ss]+vz[ss]*vz[ss];
	b2=bx2+by[ss]*by[ss]+bz[ss]*bz[ss];
	pr[ss]=(gamma-1)*(en[ss]-0.5*(ro[ss]*v2+b2));
      }

      /* Primitive variable at cell face */
#ifdef _OPENMP
#pragma omp for
#endif
      for (i=2;i<nx-2;i++){
	int sm2,sm1,ss0,sp1,sp2;
	sm2=i-2;
	sm1=i-1;
	ss0=i;
	sp1=i+1;
	sp2=i+2;
	/* Primitive variables in the stencil */
	double ros[]={ro[sm2],ro[sm1],ro[ss0],ro[sp1],ro[sp2]};
	double vxs[]={vx[sm2],vx[sm1],vx[ss0],vx[sp1],vx[sp2]};
	double vys[]={vy[sm2],vy[sm1],vy[ss0],vy[sp1],vy[sp2]};
	double vzs[]={vz[sm2],vz[sm1],vz[ss0],vz[sp1],vz[sp2]};
	double bys[]={by[sm2],by[sm1],by[ss0],by[sp1],by[sp2]};
	double bzs[]={bz[sm2],bz[sm1],bz[ss0],bz[sp1],bz[sp2]};
	double prs[]={pr[sm2],pr[sm1],pr[ss0],pr[sp1],pr[sp2]};

	double vl[7],vr[7];

	/* Characteristic decomposition for supersonic stencil */
	double cs2[]={gamma*prs[ns/2-1]/ros[ns/2-1],
		      gamma*prs[ns/2+0]/ros[ns/2+0],
		      gamma*prs[ns/2+1]/ros[ns/2+1]};
	double vv2[]={vxs[ns/2-1]*vxs[ns/2-1]+vys[ns/2-1]*vys[ns/2-1]+vzs[ns/2-1]*vzs[ns/2-1],
		      vxs[ns/2+0]*vxs[ns/2+0]+vys[ns/2+0]*vys[ns/2+0]+vzs[ns/2+0]*vzs[ns/2+0],
		      vxs[ns/2+1]*vxs[ns/2+1]+vys[ns/2+1]*vys[ns/2+1]+vzs[ns/2+1]*vzs[ns/2+1]};
	if (max(vv2[0]/cs2[0],max(vv2[1]/cs2[1],vv2[2]/cs2[2])) >= 1.0){
	  mhd_c_reconst(ros,vxs,vys,vzs,bys,bzs,prs,bx,gamma,ns,vl,vr,func_lr);
	} else{
	  mhd_a_reconst(ros,vxs,vys,vzs,bys,bzs,prs,bx,gamma,ns,vl,vr,func_lr);
	}
	
	/* Left-face @ i+1/2 */
	ss=i+1;
	ul[7*ss+0]=vl[0];	/* ro */
	ul[7*ss+1]=vl[1];	/* vx */
	ul[7*ss+2]=vl[2];	/* vy */
	ul[7*ss+3]=vl[3];	/* vz */
	ul[7*ss+4]=vl[4];	/* by */
	ul[7*ss+5]=vl[5];	/* bz */
	ul[7*ss+6]=vl[6];	/* pr */
	/* Right-face @ i-1/2 */
	ss=i;
	ur[7*ss+0]=vr[0];	/* ro */
	ur[7*ss+1]=vr[1];	/* vx */
	ur[7*ss+2]=vr[2];	/* vy */
	ur[7*ss+3]=vr[3];	/* vz */
	ur[7*ss+4]=vr[4];	/* by */
	ur[7*ss+5]=vr[5];	/* bz */
	ur[7*ss+6]=vr[6];	/* pr */
      }

      /* Numerical flux at cell face */
#ifdef _OPENMP
#pragma omp for
#endif
      for (i=3;i<nx-2;i++){
	ss=i;	
	double flux[7];
	double dvsd[2]={(vx[i]-vx[i-1]),0.0};
	func_flux(ul[7*ss+0],ul[7*ss+1],ul[7*ss+2],ul[7*ss+3],ul[7*ss+4],ul[7*ss+5],ul[7*ss+6],
		  ur[7*ss+0],ur[7*ss+1],ur[7*ss+2],ur[7*ss+3],ur[7*ss+4],ur[7*ss+5],ur[7*ss+6],
		  bx,gamma,dvsd,
		  &flux[0],&flux[1],&flux[2],&flux[3],&flux[4],&flux[5],&flux[6]);
	fx[7*ss+0]=flux[0];	/* ro */
	fx[7*ss+1]=flux[1];	/* mx */
	fx[7*ss+2]=flux[2];	/* my */
	fx[7*ss+3]=flux[3];	/* mz */
	fx[7*ss+4]=flux[4];	/* by */
	fx[7*ss+5]=flux[5];	/* bz */
	fx[7*ss+6]=flux[6];	/* en */
      }

      /* Update */
#ifdef _OPENMP
#pragma omp for
#endif
#pragma simd
      for (i=4;i<nx-4;i++){
	int sm1,ss0,sp1,sp2;
	ss=i;
	sm1=i-1;
	ss0=i;
	sp1=i+1;
	sp2=i+2;
	double du[7];
	double valx[7][4]={{fx[7*sm1+0],fx[7*ss0+0],fx[7*sp1+0],fx[7*sp2+0]},
			   {fx[7*sm1+1],fx[7*ss0+1],fx[7*sp1+1],fx[7*sp2+1]},
			   {fx[7*sm1+2],fx[7*ss0+2],fx[7*sp1+2],fx[7*sp2+2]},
			   {fx[7*sm1+3],fx[7*ss0+3],fx[7*sp1+3],fx[7*sp2+3]},
			   {fx[7*sm1+4],fx[7*ss0+4],fx[7*sp1+4],fx[7*sp2+4]},
			   {fx[7*sm1+5],fx[7*ss0+5],fx[7*sp1+5],fx[7*sp2+5]},
			   {fx[7*sm1+6],fx[7*ss0+6],fx[7*sp1+6],fx[7*sp2+6]}};
	du[0]=func_df(&valx[0][0]);
	du[1]=func_df(&valx[1][0]);
	du[2]=func_df(&valx[2][0]);
	du[3]=func_df(&valx[3][0]);
	du[4]=func_df(&valx[4][0]);
	du[5]=func_df(&valx[5][0]);
	du[6]=func_df(&valx[6][0]);
			  
	ro[ss]=rk_fac[rk][0]*ut[7*ss+0]+rk_fac[rk][1]*(ro[ss]-dtdx*du[0]);
	mx[ss]=rk_fac[rk][0]*ut[7*ss+1]+rk_fac[rk][1]*(mx[ss]-dtdx*du[1]);
	my[ss]=rk_fac[rk][0]*ut[7*ss+2]+rk_fac[rk][1]*(my[ss]-dtdx*du[2]);
	mz[ss]=rk_fac[rk][0]*ut[7*ss+3]+rk_fac[rk][1]*(mz[ss]-dtdx*du[3]);
	by[ss]=rk_fac[rk][0]*ut[7*ss+4]+rk_fac[rk][1]*(by[ss]-dtdx*du[4]);
	bz[ss]=rk_fac[rk][0]*ut[7*ss+5]+rk_fac[rk][1]*(bz[ss]-dtdx*du[5]);
	en[ss]=rk_fac[rk][0]*ut[7*ss+6]+rk_fac[rk][1]*(en[ss]-dtdx*du[6]);
      }

    } /* OpenMP */

    /* Boundary condition */
    openbc1d_mhd(ro,mx,my,mz,en,by,bz,bx,nx,xoff,gamma,0,1);
  }

  free(ut);
  free(ul);
  free(ur);
  free(fx);
  free(vx);
  free(vy);
  free(vz);
  free(pr);
}
