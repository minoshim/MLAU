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

#define ODR (2)			/* Spatial order (1,2,3,4) */
#define R_K (3)			/* Runge-Kutta order (1,2,3). Never set >3 */
#define CTW (1)			/* Flag for CT 2D upwind weighting (Minoshima+19, ApJS,242,14) */
#define HDF (1)			/* Flag for hyper diffusion (used for Hall solver) */

void mhd_fd_ct_2d(double *ro, double *mx, double *my, double *mz,
		  double *en, double *bx, double *by, double *bz,
		  double dt, double dx, double dy,
		  int nx, int ny, int xoff, int yoff, double gamma,
		  int mpi_rank, int mpi_numx, int mpi_numy);
void diff_ctfield_e(double *en, double *bx, double *by, double *bz,
		    double eta0,
		    double dt, double dx, double dy,
		    int nx, int ny, int xoff, int yoff, double gamma,
		    int mpi_rank, int mpi_numx, int mpi_numy);
void ns_viscous_2d(const double *ro, double *mx, double *my, double *mz,
		   double *en,
		   double nu0,
		   double dt, double dx, double dy,
		   int nx, int ny, int xoff, int yoff, double gamma,
		   int mpi_rank, int mpi_numx, int mpi_numy);
void hall_fd_ct_2d(double *ro, double *mx, double *my, double *mz,
		   double *en, double *bx, double *by, double *bz,
		   double dt, double dx, double dy,
		   int nx, int ny, int xoff, int yoff, double eta_h,
		   int mpi_rank, int mpi_numx, int mpi_numy);

void mpi_sdrv2d(double *f[], int nn, int nx, int ny, int xoff, int yoff, 
		int mpi_rank, int mpi_numx, int mpi_numy);
void mpi_xbc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy);
void mpi_ybc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy);

inline double calc_bc(const double *f)
/* Interpolate B-field to cell center */
{
#if (ODR <= 2)
  return(0.5*(f[1]+f[2]));
#else
  return(0.0625*(9.0*(f[1]+f[2])-(f[0]+f[3])));
#endif
}

inline void lcal_flr(const double *f, double *fl, double *fr)
/* Linear interpolation to cell face */
{
#if (ODR == 1)
  *fl=f[2]; *fr=f[2];
#elif (ODR == 2)
  /* (*fl)=0.25*(-f[1]+4.0*f[2]+f[3]); */
  /* (*fr)=0.25*(-f[3]+4.0*f[2]+f[1]); */
  (*fl)=(-f[1]+5.0*f[2]+2.0*f[3])/6.0;
  (*fr)=(-f[3]+5.0*f[2]+2.0*f[1])/6.0;
#elif (ODR == 3)
  (*fl)=0.125*(-f[1]+6.0*f[2]+3.0*f[3]);
  (*fr)=0.125*(-f[3]+6.0*f[2]+3.0*f[1]);
#else
  (*fl)=0.003125*(9.0*f[0]-56.0*f[1]+234.0*f[2]+144.0*f[3]-11.0*f[4]);
  (*fr)=0.003125*(9.0*f[4]-56.0*f[3]+234.0*f[2]+144.0*f[1]-11.0*f[0]);
#endif
}

inline void calc_flr(const double *f, double *fl, double *fr)
/* Interpolation to cell face */
{
#if (ODR == 1)
  *fl=f[2]; *fr=f[2];
#elif (ODR == 2)
  /* muscl_mm_cal_flr(f[1],f[2],f[3],fl,fr); */
  /* muscl_mc_cal_flr(f[1],f[2],f[3],fl,fr); */
  muscl_kr_cal_flr(f[1],f[2],f[3],fl,fr);
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
inline double calc_dfc(const double *f)
/* Evaluate 1st difference from co-located points */
{
#if (ODR <= 2)
  return(0.5*(f[3]-f[1]));
#else
  return((8.0*(f[3]-f[1])-(f[4]-f[0]))/12.0);
#endif
}

inline double calc_d2fc(const double *f)
/* Evaluate 2nd difference from co-located points */
{
#if (ODR <= 2)
  return((2.0*f[2]-(f[3]+f[1])));
#else
  return((-30.0*f[2]+16.0*(f[3]+f[1])-(f[4]+f[0]))/12.0);
#endif
}

void mhd_fd_ct_2d(double *ro, double *mx, double *my, double *mz,
		  double *en, double *bx, double *by, double *bz,
		  double dt, double dx, double dy,
		  int nx, int ny, int xoff, int yoff, double gamma,
		  int mpi_rank, int mpi_numx, int mpi_numy)

/* 2D finite-difference code for MHD */
/* 2nd, 3rd, and 4th order with characteristic decomposition */
/* MLAU/HLLD/ROE + Central-Upwind-CT */
/* Need four offsets (maximum) for boundary */

/* ro: density */
/* m?: moment */
/* en: energy normalized by B^2/(4 pi) */
/* b?: magnetic field */
/* gamma: specific heat ratio */

{
  int i,j,ss,rk;
  const int ns=5;		/* number of grids in stencil */
  const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const int nxy=nx*ny;
  const double dtdx=dt/dx,dtdy=dt/dy;
  double *ut,*ul,*ur,*ql,*qr,*fx,*fy;
  double *vx,*vy,*vz,*pr,*cx,*cy,*ez;
  double *ctx,*fc,eps=1e-6;
  double *m2,*dvx,*dvy;
  double (*func_bc)(const double*)=&calc_bc;
  void (*func_lr)(const double*, double*, double*)=&calc_flr;
  void (*lfun_lr)(const double*, double*, double*)=&lcal_flr;
  double (*func_df)(const double*)=&calc_df;
  void (*func_flux)(double, double, double, double, double, double, double,
		    double, double, double, double, double, double, double,
		    double, double, double*,
		    double*, double*, double*, double*, double*, double*, double*);
  func_flux=&calc_flux_mlau;
  /* func_flux=&calc_flux_hlld; */
  /* func_flux=&calc_flux_roe; */
  double *p[8];
  
  ut=(double*)malloc(sizeof(double)*8*nxy);
  ul=(double*)malloc(sizeof(double)*8*nxy);
  ur=(double*)malloc(sizeof(double)*8*nxy);
  ql=(double*)malloc(sizeof(double)*nxy);
  qr=(double*)malloc(sizeof(double)*nxy);
  fx=(double*)malloc(sizeof(double)*8*nxy);
  fy=(double*)malloc(sizeof(double)*8*nxy);
  vx=(double*)malloc(sizeof(double)*nxy);
  vy=(double*)malloc(sizeof(double)*nxy);
  vz=(double*)malloc(sizeof(double)*nxy);
  pr=(double*)malloc(sizeof(double)*nxy);
  cx=(double*)malloc(sizeof(double)*nxy);
  cy=(double*)malloc(sizeof(double)*nxy);
  ez=(double*)malloc(sizeof(double)*nxy);
  ctx=(double*)malloc(sizeof(double)*nxy);
  fc=(double*)malloc(sizeof(double)*nxy);
  m2=(double*)malloc(sizeof(double)*nxy);
  dvx=(double*)malloc(sizeof(double)*nxy);
  dvy=(double*)malloc(sizeof(double)*nxy);

  /* Copy */
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      ss=nx*j+i;
      ut[8*ss+0]=ro[ss];
      ut[8*ss+1]=mx[ss];
      ut[8*ss+2]=my[ss];
      ut[8*ss+3]=mz[ss];
      ut[8*ss+4]=bx[ss];
      ut[8*ss+5]=by[ss];
      ut[8*ss+6]=bz[ss];
      ut[8*ss+7]=en[ss];
    }
  }

  for (rk=0;rk<R_K;rk++){
    
#ifdef _OPENMP
#pragma omp parallel private(i,j,ss)
#endif
    {

      /* Calculate B at cell center */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=0;j<ny;j++){
#pragma simd
	for (i=1;i<nx-2;i++){
	  double val[]={bx[nx*j+(i-1)],bx[nx*j+(i+0)],
			bx[nx*j+(i+1)],bx[nx*j+(i+2)]};
	  cx[nx*j+i]=func_bc(val);
	}
      }
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=1;j<ny-2;j++){
#pragma simd
	for (i=0;i<nx;i++){
	  double val[]={by[nx*(j-1)+i],by[nx*(j+0)+i],
			by[nx*(j+1)+i],by[nx*(j+2)+i]};
	  cy[nx*j+i]=func_bc(val);
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	p[0]=cx;
	p[1]=cy;
	mpi_sdrv2d(p,2,nx,ny,xoff,yoff,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(&cx[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(&cy[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(&cx[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(&cy[0],nx,ny,xoff,yoff,0,-1,mpi_rank,mpi_numx,mpi_numy);
      }

      /* Primitive variable at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
	double iro,v2,b2;
#pragma simd
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  iro=1.0/ro[ss];
	  vx[ss]=mx[ss]*iro;
	  vy[ss]=my[ss]*iro;
	  vz[ss]=mz[ss]*iro;
	  v2=(vx[ss]*vx[ss]+vy[ss]*vy[ss]+vz[ss]*vz[ss]);
	  b2=(cx[ss]*cx[ss]+cy[ss]*cy[ss]+bz[ss]*bz[ss]);
	  pr[ss]=(gamma-1)*(en[ss]-0.5*(ro[ss]*v2+b2));
	  m2[ss]=v2/(gamma*pr[ss]*iro);
	  ez[ss]=0.0;		/* Necessary initialize at cell corner */
	  ctx[ss]=0.5;
	}
      }
      /* CT E-field 2D upwind weighting */
#if (CTW)
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=1;j<ny;j++){
	double vxc,vyc,roc,bxc,byc;
#pragma simd
	for (i=1;i<nx;i++){
	  /* Upwinding with respect to Alfven mode */
	  vxc=0.25*(vx[nx*(j-1)+(i-1)]+vx[nx*(j-1)+i]+
		    vx[nx*j+(i-1)]+vx[nx*j+i]);
	  vyc=0.25*(vy[nx*(j-1)+(i-1)]+vy[nx*(j-1)+i]+
		    vy[nx*j+(i-1)]+vy[nx*j+i]);
	  roc=0.25*(ro[nx*(j-1)+(i-1)]+ro[nx*(j-1)+i]+
		    ro[nx*j+(i-1)]+ro[nx*j+i]);
	  roc=1.0/sqrt(roc);
	  bxc=0.5*(bx[nx*(j-1)+i]+bx[nx*j+i]);
	  byc=0.5*(by[nx*j+(i-1)]+by[nx*j+i]);
	  vxc=fabs(vxc)+fabs(bxc*roc);
	  vyc=fabs(vyc)+fabs(byc*roc);

	  ctx[nx*j+i]=(vxc+0.5*eps)/(vxc+vyc+eps);
	}
      }
#endif

      /* Calculate dvx and dvy for shock detection */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=1;j<ny-1;j++){
#pragma simd
	for (i=1;i<nx-1;i++){
	  ss=nx*j+i;
	  dvx[ss]=min((vx[nx*j+(i+1)]-vx[nx*j+i]),
		      (vx[nx*j+i]-vx[nx*j+(i-1)]));
	  dvy[ss]=min((vy[nx*(j+1)+i]-vy[nx*j+i]),
		      (vy[nx*j+i]-vy[nx*(j-1)+i]));
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	p[0]=dvx;
	p[1]=dvy;
	mpi_sdrv2d(p,2,nx,ny,xoff,yoff,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(&dvx[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(&dvy[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(&dvx[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(&dvy[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
      }

      /* Primitive variable at cell face along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
	int sl,sr;
	int sm2,sm,sp,sp2;
	for (i=2;i<nx-2;i++){
	  ss=nx*j+i;
	  sl=nx*j+(i+1);
	  sr=nx*j+i;
	  sm2=nx*j+(i-2);
	  sm=nx*j+(i-1);
	  sp=nx*j+(i+1);
	  sp2=nx*j+(i+2);
	  /* Primitive variables in the stencil */
	  double ros[]={ro[sm2],ro[sm],ro[ss],ro[sp],ro[sp2]};
	  double vxs[]={vx[sm2],vx[sm],vx[ss],vx[sp],vx[sp2]};
	  double vys[]={vy[sm2],vy[sm],vy[ss],vy[sp],vy[sp2]};
	  double vzs[]={vz[sm2],vz[sm],vz[ss],vz[sp],vz[sp2]};
	  double bxs[]={cx[sm2],cx[sm],cx[ss],cx[sp],cx[sp2]};
	  double bys[]={cy[sm2],cy[sm],cy[ss],cy[sp],cy[sp2]};
	  double bzs[]={bz[sm2],bz[sm],bz[ss],bz[sp],bz[sp2]};
	  double prs[]={pr[sm2],pr[sm],pr[ss],pr[sp],pr[sp2]};

	  double vl[7],vr[7];

	  if (max(m2[sm],max(m2[ss],m2[sp])) >= 1.0){
	    mhd_c_reconst(ros,vxs,vys,vzs,bys,bzs,prs,bxs[ns/2],gamma,ns,vl,vr,func_lr);
	  } else{
	    mhd_a_reconst(ros,vxs,vys,vzs,bys,bzs,prs,bxs[ns/2],gamma,ns,vl,vr,func_lr);
	  }

	  /* Left-face @ i+1/2 */
	  ul[8*sl+0]=vl[0];	/* ro */
	  ul[8*sl+1]=vl[1];	/* vx */
	  ul[8*sl+2]=vl[2];	/* vy */
	  ul[8*sl+3]=vl[3];	/* vz */
	  ul[8*sl+4]=bx[sl];	/* bx */
	  ul[8*sl+5]=vl[4];	/* by */
	  ul[8*sl+6]=vl[5];	/* bz */
	  ul[8*sl+7]=vl[6];	/* pr */
	  /* Right-face @ i-1/2 */
	  ur[8*sr+0]=vr[0];	/* ro */
	  ur[8*sr+1]=vr[1];	/* vx */
	  ur[8*sr+2]=vr[2];	/* vy */
	  ur[8*sr+3]=vr[3];	/* vz */
	  ur[8*sr+4]=bx[sr];	/* bx */
	  ur[8*sr+5]=vr[4];	/* by */
	  ur[8*sr+6]=vr[5];	/* bz */
	  ur[8*sr+7]=vr[6];	/* pr */

	  /* Linear interpolation of numerical flux of By */
	  double val[]={cy[sm2]*vx[sm2],cy[sm]*vx[sm],cy[ss]*vx[ss],cy[sp]*vx[sp],cy[sp2]*vx[sp2]};
	  lfun_lr(val,&vl[0],&vr[0]);
	  lfun_lr(vys,&vl[1],&vr[1]);
	  ql[sl]=vl[0]-bx[sl]*vl[1];
	  qr[sr]=vr[0]-bx[sr]*vr[1];
	}
      }
      /* Numerical flux at cell face along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
#pragma simd
	for (i=3;i<nx-2;i++){
	  double bn,flux[8];
	  ss=nx*j+i;
	  bn=0.5*(ul[8*ss+4]+ur[8*ss+4]);

	  double dvsd[2]={(vx[nx*j+i]-vx[nx*j+(i-1)]),min(dvy[nx*j+(i-1)],dvy[nx*j+i])};
	  func_flux(ul[8*ss+0],ul[8*ss+1],ul[8*ss+2],ul[8*ss+3],ul[8*ss+5],ul[8*ss+6],ul[8*ss+7],
		    ur[8*ss+0],ur[8*ss+1],ur[8*ss+2],ur[8*ss+3],ur[8*ss+5],ur[8*ss+6],ur[8*ss+7],
		    bn,gamma,dvsd,
		    &flux[0],&flux[1],&flux[2],&flux[3],&flux[5],&flux[6],&flux[7]);
	  
	  fx[8*ss+0]=flux[0];	/* ro */
	  fx[8*ss+1]=flux[1];	/* mx */
	  fx[8*ss+2]=flux[2];	/* my */
	  fx[8*ss+3]=flux[3];	/* mz */
	  fx[8*ss+4]=0;		/* bx */
	  fx[8*ss+5]=flux[5];	/* by */
	  fx[8*ss+6]=flux[6];	/* bz */
	  fx[8*ss+7]=flux[7];	/* en */
	  /* Divide central and upwind parts in numerical flux of By */
	  fc[ss]=0.5*(ql[ss]+qr[ss]); /* Central part */
	  fx[8*ss+5]-=fc[ss];	      /* Upwind part */
	}
      }

      /* Numerical flux of By at cell corner along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-2;j++){
	int sl,sr;
	int sm2,sm,sp,sp2;
#pragma simd
	for (i=3;i<nx-2;i++){
	  ss=nx*j+i;
	  sl=nx*(j+1)+i;
	  sr=nx*j+i;
	  sm2=nx*(j-2)+i;
	  sm=nx*(j-1)+i;
	  sp=nx*(j+1)+i;
	  sp2=nx*(j+2)+i;
	  double val[]={fx[8*sm2+5],fx[8*sm+5],fx[8*ss+5],fx[8*sp+5],fx[8*sp2+5]};
	  double vac[]={fc[sm2],fc[sm],fc[ss],fc[sp],fc[sp2]};

	  double vl[2],vr[2];
	  lfun_lr(val,&vl[0],&vr[0]); /* Upwind part */
	  lfun_lr(vac,&vl[1],&vr[1]); /* Central part */
	  ul[2*sl+0]=vl[0];
	  ur[2*sr+0]=vr[0];
	  ul[2*sl+1]=vl[1];
	  ur[2*sr+1]=vr[1];
	}
      }
      /* E-field at cell corner along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
#pragma simd
	for (i=3;i<nx-2;i++){
	  ss=nx*j+i;
	  ez[ss]+=-0.5*((ul[2*ss+0]+ur[2*ss+0])+(1.0-ctx[ss])*(ul[2*ss+1]+ur[2*ss+1]));
	}
      }

      /* Primitive variable at cell face along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-2;j++){
	int sl,sr;
	int sm2,sm,sp,sp2;
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  sl=nx*(j+1)+i;
	  sr=nx*j+i;
	  sm2=nx*(j-2)+i;
	  sm=nx*(j-1)+i;
	  sp=nx*(j+1)+i;
	  sp2=nx*(j+2)+i;
	  /* Primitive variables in the stencil */
	  double ros[]={ro[sm2],ro[sm],ro[ss],ro[sp],ro[sp2]};
	  double vxs[]={vx[sm2],vx[sm],vx[ss],vx[sp],vx[sp2]};
	  double vys[]={vy[sm2],vy[sm],vy[ss],vy[sp],vy[sp2]};
	  double vzs[]={vz[sm2],vz[sm],vz[ss],vz[sp],vz[sp2]};
	  double bxs[]={cx[sm2],cx[sm],cx[ss],cx[sp],cx[sp2]};
	  double bys[]={cy[sm2],cy[sm],cy[ss],cy[sp],cy[sp2]};
	  double bzs[]={bz[sm2],bz[sm],bz[ss],bz[sp],bz[sp2]};
	  double prs[]={pr[sm2],pr[sm],pr[ss],pr[sp],pr[sp2]};

	  double vl[7],vr[7];

	  if (max(m2[sm],max(m2[ss],m2[sp])) >= 1.0){
	    mhd_c_reconst(ros,vys,vzs,vxs,bzs,bxs,prs,bys[ns/2],gamma,ns,vl,vr,func_lr);
	  } else{
	    mhd_a_reconst(ros,vys,vzs,vxs,bzs,bxs,prs,bys[ns/2],gamma,ns,vl,vr,func_lr);
	  }

	  /* Left-face @ j+1/2 */
	  ul[8*sl+0]=vl[0];	/* ro */
	  ul[8*sl+1]=vl[1];	/* vy */
	  ul[8*sl+2]=vl[2];	/* vz */
	  ul[8*sl+3]=vl[3];	/* vx */
	  ul[8*sl+4]=by[sl];	/* by */
	  ul[8*sl+5]=vl[4];	/* bz */
	  ul[8*sl+6]=vl[5];	/* bx */
	  ul[8*sl+7]=vl[6];	/* pr */
	  /* Right-face @ j-1/2 */
	  ur[8*sr+0]=vr[0];	/* ro */
	  ur[8*sr+1]=vr[1];	/* vy */
	  ur[8*sr+2]=vr[2];	/* vz */
	  ur[8*sr+3]=vr[3];	/* vx */
	  ur[8*sr+4]=by[sr];	/* by */
	  ur[8*sr+5]=vr[4];	/* bz */
	  ur[8*sr+6]=vr[5];	/* bx */
	  ur[8*sr+7]=vr[6];	/* pr */

	  /* Linear interpolation of numerical flux of Bx */
	  double val[]={cx[sm2]*vy[sm2],cx[sm]*vy[sm],cx[ss]*vy[ss],cx[sp]*vy[sp],cx[sp2]*vy[sp2]};
	  lfun_lr(val,&vl[0],&vr[0]);
	  lfun_lr(vxs,&vl[1],&vr[1]);
	  ql[sl]=vl[0]-by[sl]*vl[1];
	  qr[sr]=vr[0]-by[sr]*vr[1];
	}
      }
      /* Numerical flux at cell face along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
#pragma simd
	for (i=0;i<nx;i++){
	  double bn,flux[8];
	  ss=nx*j+i;
	  bn=0.5*(ul[8*ss+4]+ur[8*ss+4]);

	  double dvsd[2]={(vy[nx*j+i]-vy[nx*(j-1)+i]),min(dvx[nx*(j-1)+i],dvx[nx*j+i])};
	  func_flux(ul[8*ss+0],ul[8*ss+1],ul[8*ss+2],ul[8*ss+3],ul[8*ss+5],ul[8*ss+6],ul[8*ss+7],
		    ur[8*ss+0],ur[8*ss+1],ur[8*ss+2],ur[8*ss+3],ur[8*ss+5],ur[8*ss+6],ur[8*ss+7],
		    bn,gamma,dvsd,
		    &flux[0],&flux[1],&flux[2],&flux[3],&flux[5],&flux[6],&flux[7]);

	  fy[8*ss+0]=flux[0];	/* ro */
	  fy[8*ss+2]=flux[1];	/* my */
	  fy[8*ss+3]=flux[2];	/* mz */
	  fy[8*ss+1]=flux[3];	/* mx */
	  fy[8*ss+5]=0;		/* by */
	  fy[8*ss+6]=flux[5];	/* bz */
	  fy[8*ss+4]=flux[6];	/* bx */
	  fy[8*ss+7]=flux[7];	/* en */
	  /* Divide central and upwind parts in numerical flux of Bx */
	  fc[ss]=0.5*(ql[ss]+qr[ss]); /* Central part */
	  fy[8*ss+4]-=fc[ss];	      /* Upwind part */
	}
      }

      /* Numerical flux of Bx at cell corner along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
	int sl,sr;
	int sm2,sm,sp,sp2;
#pragma simd
	for (i=2;i<nx-2;i++){
	  ss=nx*j+i;
	  sl=nx*j+(i+1);
	  sr=nx*j+i;
	  sm2=nx*j+(i-2);
	  sm=nx*j+(i-1);
	  sp=nx*j+(i+1);
	  sp2=nx*j+(i+2);
	  double val[]={fy[8*sm2+4],fy[8*sm+4],fy[8*ss+4],fy[8*sp+4],fy[8*sp2+4]};
	  double vac[]={fc[sm2],fc[sm],fc[ss],fc[sp],fc[sp2]};

	  double vl[2],vr[2];
	  lfun_lr(val,&vl[0],&vr[0]); /* Upwind part */
	  lfun_lr(vac,&vl[1],&vr[1]); /* Central part */
	  ul[2*sl+0]=vl[0];
	  ur[2*sr+0]=vr[0];
	  ul[2*sl+1]=vl[1];
	  ur[2*sr+1]=vr[1];
	}
      }
      /* E-field at cell corner along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
#pragma simd
	for (i=3;i<nx-2;i++){
	  ss=nx*j+i;
	  ez[ss]+=+0.5*((ul[2*ss+0]+ur[2*ss+0])+ctx[ss]*(ul[2*ss+1]+ur[2*ss+1]));
	}
      }

      /* Update variable at cell center */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=yoff;j<ny-yoff;j++){
	int sip,sjp,si0,sj0;
	int sim,sip2,sjm,sjp2;
#pragma simd
	for (i=xoff;i<nx-xoff;i++){
	  ss=nx*j+i;
	  sip=nx*j+(i+1);
	  sjp=nx*(j+1)+i;
	  si0=nx*j+i;
	  sj0=nx*j+i;
	  sim=nx*j+(i-1);
	  sjm=nx*(j-1)+i;
	  sip2=nx*j+(i+2);
	  sjp2=nx*(j+2)+i;
	  double du[2][8];
	  /* dF/dx */
	  double valx[8][4]={{fx[8*sim+0],fx[8*si0+0],fx[8*sip+0],fx[8*sip2+0]},
			     {fx[8*sim+1],fx[8*si0+1],fx[8*sip+1],fx[8*sip2+1]},
			     {fx[8*sim+2],fx[8*si0+2],fx[8*sip+2],fx[8*sip2+2]},
			     {fx[8*sim+3],fx[8*si0+3],fx[8*sip+3],fx[8*sip2+3]},
			     {fx[8*sim+4],fx[8*si0+4],fx[8*sip+4],fx[8*sip2+4]},
			     {fx[8*sim+5],fx[8*si0+5],fx[8*sip+5],fx[8*sip2+5]},
			     {fx[8*sim+6],fx[8*si0+6],fx[8*sip+6],fx[8*sip2+6]},
			     {fx[8*sim+7],fx[8*si0+7],fx[8*sip+7],fx[8*sip2+7]}};
	  du[0][0]=func_df(&valx[0][0]);
	  du[0][1]=func_df(&valx[1][0]);
	  du[0][2]=func_df(&valx[2][0]);
	  du[0][3]=func_df(&valx[3][0]);
	  du[0][6]=func_df(&valx[6][0]);
	  du[0][7]=func_df(&valx[7][0]);
	  /* dG/dy */
	  double valy[8][4]={{fy[8*sjm+0],fy[8*sj0+0],fy[8*sjp+0],fy[8*sjp2+0]},
			     {fy[8*sjm+1],fy[8*sj0+1],fy[8*sjp+1],fy[8*sjp2+1]},
			     {fy[8*sjm+2],fy[8*sj0+2],fy[8*sjp+2],fy[8*sjp2+2]},
			     {fy[8*sjm+3],fy[8*sj0+3],fy[8*sjp+3],fy[8*sjp2+3]},
			     {fy[8*sjm+4],fy[8*sj0+4],fy[8*sjp+4],fy[8*sjp2+4]},
			     {fy[8*sjm+5],fy[8*sj0+5],fy[8*sjp+5],fy[8*sjp2+5]},
			     {fy[8*sjm+6],fy[8*sj0+6],fy[8*sjp+6],fy[8*sjp2+6]},
			     {fy[8*sjm+7],fy[8*sj0+7],fy[8*sjp+7],fy[8*sjp2+7]}};
	  du[1][0]=func_df(&valy[0][0]);
	  du[1][1]=func_df(&valy[1][0]);
	  du[1][2]=func_df(&valy[2][0]);
	  du[1][3]=func_df(&valy[3][0]);
	  du[1][6]=func_df(&valy[6][0]);
	  du[1][7]=func_df(&valy[7][0]);

	  ro[ss]=rk_fac[rk][0]*ut[8*ss+0]+rk_fac[rk][1]*(ro[ss]-dtdx*du[0][0]-dtdy*du[1][0]);
	  mx[ss]=rk_fac[rk][0]*ut[8*ss+1]+rk_fac[rk][1]*(mx[ss]-dtdx*du[0][1]-dtdy*du[1][1]);
	  my[ss]=rk_fac[rk][0]*ut[8*ss+2]+rk_fac[rk][1]*(my[ss]-dtdx*du[0][2]-dtdy*du[1][2]);
	  mz[ss]=rk_fac[rk][0]*ut[8*ss+3]+rk_fac[rk][1]*(mz[ss]-dtdx*du[0][3]-dtdy*du[1][3]);
	  bz[ss]=rk_fac[rk][0]*ut[8*ss+6]+rk_fac[rk][1]*(bz[ss]-dtdx*du[0][6]-dtdy*du[1][6]);
	  en[ss]=rk_fac[rk][0]*ut[8*ss+7]+rk_fac[rk][1]*(en[ss]-dtdx*du[0][7]-dtdy*du[1][7]);
	}
      }
      /* Update CT Bx */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=yoff;j<ny-yoff;j++){
	int sp,s0;
	int sm,sp2;
	double du;
#pragma simd
	for (i=xoff;i<nx-xoff+1;i++){
	  ss=nx*j+i;
	  sp=nx*(j+1)+i;
	  s0=nx*j+i;
	  sm=nx*(j-1)+i;
	  sp2=nx*(j+2)+i;
	  double val[]={ez[sm],ez[s0],ez[sp],ez[sp2]};
	  du=func_df(val);
	  bx[ss]=rk_fac[rk][0]*ut[8*ss+4]+rk_fac[rk][1]*(bx[ss]-dtdy*du);
	}
      }
      /* Update CT By */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=yoff;j<ny-yoff+1;j++){
	int sp,s0;
	int sm,sp2;
	double du;
#pragma simd
	for (i=xoff;i<nx-xoff;i++){
	  ss=nx*j+i;
	  sp=nx*j+(i+1);
	  s0=nx*j+i;
	  sm=nx*j+(i-1);
	  sp2=nx*j+(i+2);
	  double val[]={ez[sm],ez[s0],ez[sp],ez[sp2]};
	  du=func_df(val);
	  by[ss]=rk_fac[rk][0]*ut[8*ss+5]+rk_fac[rk][1]*(by[ss]+dtdx*du);
	}
      }

    } /* OpenMP */

    /* Boundary condition */
    p[0]=ro;
    p[1]=mx;
    p[2]=my;
    p[3]=mz;
    p[4]=bx;
    p[5]=by;
    p[6]=bz;
    p[7]=en;
    mpi_sdrv2d(p,8,nx,ny,xoff,yoff,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(&ro[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(&mx[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(&my[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(&mz[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(&en[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(&bx[0],nx,ny,xoff,yoff,1,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(&by[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(&bz[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&ro[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&mx[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&my[0],nx,ny,xoff,yoff,0,-1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&mz[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&en[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&bx[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&by[0],nx,ny,xoff,yoff,1,-1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&bz[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
  }

  free(ut);
  free(ul);
  free(ur);
  free(ql);
  free(qr);
  free(fx);
  free(fy);
  free(vx);
  free(vy);
  free(vz);
  free(pr);
  free(cx);
  free(cy);
  free(ez);
  free(ctx);
  free(fc);
  free(m2);
  free(dvx);
  free(dvy);
}

void diff_ctfield_e(double *en, double *bx, double *by, double *bz,
		    double eta0,
		    double dt, double dx, double dy,
		    int nx, int ny, int xoff, int yoff, double gamma,
		    int mpi_rank, int mpi_numx, int mpi_numy)
/* Magnetic field diffusion with uniform resistivity */
/* Use CT grid spacing */
/* Update energy as well as mag. field */
{
  int i,j,ss,rk;
  const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const int nxy=nx*ny;
  const double dtdx=dt/dx,dtdy=dt/dy;
  const double idx=1./dx,idy=1./dy;
  double *ut,*fx,*fy;
  double *ex,*ey,*ez;
  double (*func_bc)(const double*)=&calc_bc;
  double (*func_df)(const double*)=&calc_df;
  double *p[3];

  ut=(double*)malloc(sizeof(double)*4*nxy);
  fx=(double*)malloc(sizeof(double)*nxy);
  fy=(double*)malloc(sizeof(double)*nxy);
  ex=(double*)malloc(sizeof(double)*nxy);
  ey=(double*)malloc(sizeof(double)*nxy);
  ez=(double*)malloc(sizeof(double)*nxy);

  /* Copy */
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      ss=nx*j+i;
      ut[4*ss+0]=bx[ss];
      ut[4*ss+1]=by[ss];
      ut[4*ss+2]=bz[ss];
      ut[4*ss+3]=en[ss];
    }
  }

  for (rk=0;rk<R_K;rk++){

#ifdef _OPENMP
#pragma omp parallel private(i,j,ss)
#endif
    {

      /* Calculate resistive e-field */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-1;j++){
	int ssi0,ssi1,ssi2,ssi3,ssj0,ssj1,ssj2,ssj3;
	double jx,jy,jz;
#pragma simd
	for (i=2;i<nx-1;i++){
	  ss=nx*j+i;
	  ssi0=nx*j+(i-2);
	  ssi1=nx*j+(i-1);
	  ssi2=nx*j+(i+0);
	  ssi3=nx*j+(i+1);
	  ssj0=nx*(j-2)+i;
	  ssj1=nx*(j-1)+i;
	  ssj2=nx*(j+0)+i;
	  ssj3=nx*(j+1)+i;

	  double bxj[]={bx[ssj0],bx[ssj1],bx[ssj2],bx[ssj3]};
	  double byi[]={by[ssi0],by[ssi1],by[ssi2],by[ssi3]};
	  double bzi[]={bz[ssi0],bz[ssi1],bz[ssi2],bz[ssi3]};
	  double bzj[]={bz[ssj0],bz[ssj1],bz[ssj2],bz[ssj3]};
	  jx=+func_df(bzj)*idy;
	  jy=-func_df(bzi)*idx;
	  jz=+func_df(byi)*idx-func_df(bxj)*idy;

	  ex[ss]=eta0*jx;
	  ey[ss]=eta0*jy;
	  ez[ss]=eta0*jz;
	}
      }
      /* Calculate numerical flux for energy along X */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=2;j<ny-2;j++){
	int ssi0,ssi1,ssi2,ssi3;
	double ep,bp,em,bm;
#pragma simd
	for (i=2;i<nx-1;i++){
	  ss=nx*j+i;
	  ssi0=nx*j+(i-2);
	  ssi1=nx*j+(i-1);
	  ssi2=nx*j+(i+0);
	  ssi3=nx*j+(i+1);

	  double byi[]={by[ssi0],by[ssi1],by[ssi2],by[ssi3]};
	  double byp[]={by[nx*(j+1)+(i-2)],by[nx*(j+1)+(i-1)],
			by[nx*(j+1)+(i+0)],by[nx*(j+1)+(i+1)]};
	  double bzi[]={bz[ssi0],bz[ssi1],bz[ssi2],bz[ssi3]};
	  ep=ey[ss];
	  bp=func_bc(bzi);
	  em=0.5*(ez[ss]+ez[nx*(j+1)+i]); /* Use 2nd order to orthogonal direc. */
	  bm=0.5*(func_bc(byi)+func_bc(byp));

	  fx[ss]=ep*bp-em*bm;
	}
      }
      /* Calculate numerical flux for energy along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-1;j++){
	int ssj0,ssj1,ssj2,ssj3;
	double ep,bp,em,bm;
#pragma simd
	for (i=2;i<nx-2;i++){
	  ss=nx*j+i;
	  ssj0=nx*(j-2)+i;
	  ssj1=nx*(j-1)+i;
	  ssj2=nx*(j+0)+i;
	  ssj3=nx*(j+1)+i;

	  double bxj[]={bx[ssj0],bx[ssj1],bx[ssj2],bx[ssj3]};
	  double bxp[]={bx[nx*(j-2)+(i+1)],bx[nx*(j-1)+(i+1)],
			bx[nx*(j+0)+(i+1)],bx[nx*(j+1)+(i+1)]};
	  double bzj[]={bz[ssj0],bz[ssj1],bz[ssj2],bz[ssj3]};
	  ep=0.5*(ez[ss]+ez[nx*j+(i+1)]); /* Use 2nd order to orthogonal direc. */
	  bp=0.5*(func_bc(bxj)+func_bc(bxp));
	  em=ex[ss];
	  bm=func_bc(bzj);

	  fy[ss]=ep*bp-em*bm;
	}
      }

      /* Update CT Bx */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=yoff;j<ny-yoff;j++){
	int sp;
	int sm,sp2;
	double du;
#pragma simd
	for (i=xoff;i<nx-xoff+1;i++){
	  ss=nx*j+i;
	  sp=nx*(j+1)+i;
	  sm=nx*(j-1)+i;
	  sp2=nx*(j+2)+i;
	  double val[]={ez[sm],ez[ss],ez[sp],ez[sp2]};
	  du=func_df(val);
	  bx[ss]=rk_fac[rk][0]*ut[4*ss+0]+rk_fac[rk][1]*(bx[ss]-dtdy*du);
	}
      }
      /* Update CT By */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=yoff;j<ny-yoff+1;j++){
	int sp;
	int sm,sp2;
	double du;
#pragma simd
	for (i=xoff;i<nx-xoff;i++){
	  ss=nx*j+i;
	  sp=nx*j+(i+1);
	  sm=nx*j+(i-1);
	  sp2=nx*j+(i+2);
	  double val[]={ez[sm],ez[ss],ez[sp],ez[sp2]};
	  du=func_df(val);
	  by[ss]=rk_fac[rk][0]*ut[4*ss+1]+rk_fac[rk][1]*(by[ss]+dtdx*du);
	}
      }
      /* Update Bz and energy @ center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=yoff;j<ny-yoff;j++){
	int sip,sjp;
	int sim,sip2,sjm,sjp2;
#pragma simd
	for (i=xoff;i<nx-xoff;i++){
	  ss=nx*j+i;
	  sip=nx*j+(i+1);
	  sjp=nx*(j+1)+i;
	  sim=nx*j+(i-1);
	  sjm=nx*(j-1)+i;
	  sip2=nx*j+(i+2);
	  sjp2=nx*(j+2)+i;
	  double du[2][2];
	  /* dF/dx */
	  double valx[2][4]={{ey[sim],ey[ss],ey[sip],ey[sip2]},
			     {fx[sim],fx[ss],fx[sip],fx[sip2]}};
	  du[0][0]=func_df(&valx[0][0]);
	  du[0][1]=func_df(&valx[1][0]);
	  /* dG/dy */
	  double valy[2][4]={{ex[sjm],ex[ss],ex[sjp],ex[sjp2]},
			     {fy[sjm],fy[ss],fy[sjp],fy[sjp2]}};
	  du[1][0]=func_df(&valy[0][0]);
	  du[1][1]=func_df(&valy[1][0]);

	  bz[ss]=rk_fac[rk][0]*ut[4*ss+2]+rk_fac[rk][1]*(bz[ss]-dtdx*du[0][0]+dtdy*du[1][0]);
	  en[ss]=rk_fac[rk][0]*ut[4*ss+3]+rk_fac[rk][1]*(en[ss]-dtdx*du[0][1]-dtdy*du[1][1]);
	}
      }
    } /* OpenMP */

    p[0]=bx;
    p[1]=by;
    p[2]=bz;
    mpi_sdrv2d(p,3,nx,ny,xoff,yoff,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(bx,nx,ny,xoff,yoff,1,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(by,nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(bz,nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(bx,nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(by,nx,ny,xoff,yoff,1,-1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(bz,nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
  }
  p[0]=en;
  mpi_sdrv2d(p,1,nx,ny,xoff,yoff,mpi_rank,mpi_numx,mpi_numy);
  mpi_xbc2d(en,nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
  mpi_ybc2d(en,nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);

  free(ut);
  free(fx);
  free(fy);
  free(ex);
  free(ey);
  free(ez);
}

void ns_viscous_2d(const double *ro, double *mx, double *my, double *mz,
		   double *en,
		   double nu0,
		   double dt, double dx, double dy,
		   int nx, int ny, int xoff, int yoff, double gamma,
		   int mpi_rank, int mpi_numx, int mpi_numy)
/* Solve viscous term of 2D Navier-Stokes equation */
/* nu0 is kinetic viscosity (assumed constant) */
/* Explicit, MPI parallel code */
{
  int i,j,ss,rk;
  const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const int nxy=nx*ny;
  const double idx=1./dx,idy=1./dy;
  const double dtdx=dt*idx,dtdy=dt*idy;
  const double c2=2./3.,c4=4./3.;
  double *ir,*ut,*fx,*fy;
  double (*func_bc)(const double*)=&calc_bc;
  double (*func_df)(const double*)=&calc_df;
  double *p[3];

  ir=(double*)malloc(sizeof(double)*nxy);
  ut=(double*)malloc(sizeof(double)*4*nxy);
  fx=(double*)malloc(sizeof(double)*4*nxy);
  fy=(double*)malloc(sizeof(double)*4*nxy);

  /* Copy */
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      ss=nx*j+i;
      ir[ss]=1.0/ro[ss];
      mx[ss]*=ir[ss];		/* mx => vx */
      my[ss]*=ir[ss];		/* my => vy */
      mz[ss]*=ir[ss];		/* mz => vz */
      ut[4*ss+0]=mx[ss];
      ut[4*ss+1]=my[ss];
      ut[4*ss+2]=mz[ss];
      ut[4*ss+3]=en[ss];
    }
  }

  for (rk=0;rk<R_K;rk++){
#ifdef _OPENMP
#pragma omp parallel private(i,j,ss)
#endif
    {

      /* Numerical flux at cell face along X */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=1;j<ny-1;j++){
	int sm2,sm,s0,sp;
	double ro0,rnu,vx0,vy0,vz0;
	double dvxdx,dvxdy,dvydx,dvydy,dvzdx;
#pragma simd
	for (i=2;i<nx-1;i++){
	  double dv[4];
	  ss=nx*j+i;
	  sm2=nx*j+(i-2);
	  sm=nx*j+(i-1);
	  s0=nx*j+(i+0);
	  sp=nx*j+(i+1);

	  ro0=0.5*(ro[s0]+ro[sm]); /* Use 2nd order to avoid negative */
	  rnu=ro0*nu0;

	  double vxtmp[4]={mx[sm2],mx[sm],mx[s0],mx[sp]};
	  double vytmp[4]={my[sm2],my[sm],my[s0],my[sp]};
	  double vztmp[4]={mz[sm2],mz[sm],mz[s0],mz[sp]};
	  vx0=func_bc(vxtmp);
	  vy0=func_bc(vytmp);
	  vz0=func_bc(vztmp);

	  dvxdx=func_df(vxtmp)*idx;
	  dvydx=func_df(vytmp)*idx;
	  dvzdx=func_df(vztmp)*idx;

	  dv[0]=0.5*(mx[nx*(j+1)+(i-2)]-mx[nx*(j-1)+(i-2)]);
	  dv[1]=0.5*(mx[nx*(j+1)+(i-1)]-mx[nx*(j-1)+(i-1)]);
	  dv[2]=0.5*(mx[nx*(j+1)+(i+0)]-mx[nx*(j-1)+(i+0)]);
	  dv[3]=0.5*(mx[nx*(j+1)+(i+1)]-mx[nx*(j-1)+(i+1)]);
	  dvxdy=func_bc(dv)*idy;

	  dv[0]=0.5*(my[nx*(j+1)+(i-2)]-my[nx*(j-1)+(i-2)]);
	  dv[1]=0.5*(my[nx*(j+1)+(i-1)]-my[nx*(j-1)+(i-1)]);
	  dv[2]=0.5*(my[nx*(j+1)+(i+0)]-my[nx*(j-1)+(i+0)]);
	  dv[3]=0.5*(my[nx*(j+1)+(i+1)]-my[nx*(j-1)+(i+1)]);
	  dvydy=func_bc(dv)*idy;

	  fx[4*ss+0]=-rnu*(c4*dvxdx-c2*dvydy);
	  fx[4*ss+1]=-rnu*(dvydx+dvxdy);
	  fx[4*ss+2]=-rnu*(dvzdx);
	  fx[4*ss+3]=fx[4*ss+0]*vx0+fx[4*ss+1]*vy0+fx[4*ss+2]*vz0;
	}
      }

      /* Numerical flux at cell face along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-1;j++){
	int sm2,sm,s0,sp;
	double ro0,rnu,vx0,vy0,vz0;
	double dvxdx,dvxdy,dvydx,dvydy,dvzdy;
#pragma simd
	for (i=1;i<nx-1;i++){
	  double dv[4];
	  ss=nx*j+i;
	  sm2=nx*(j-2)+i;
	  sm=nx*(j-1)+i;
	  s0=nx*(j+0)+i;
	  sp=nx*(j+1)+i;
	  
	  ro0=0.5*(ro[s0]+ro[sm]); /* Use 2nd order to avoid negative */
	  rnu=ro0*nu0;

	  double vxtmp[4]={mx[sm2],mx[sm],mx[s0],mx[sp]};
	  double vytmp[4]={my[sm2],my[sm],my[s0],my[sp]};
	  double vztmp[4]={mz[sm2],mz[sm],mz[s0],mz[sp]};
	  vx0=func_bc(vxtmp);
	  vy0=func_bc(vytmp);
	  vz0=func_bc(vztmp);

	  dvxdy=func_df(vxtmp)*idy;
	  dvydy=func_df(vytmp)*idy;
	  dvzdy=func_df(vztmp)*idy;

	  dv[0]=0.5*(mx[nx*(j-2)+(i+1)]-mx[nx*(j-2)+(i-1)]);
	  dv[1]=0.5*(mx[nx*(j-1)+(i+1)]-mx[nx*(j-1)+(i-1)]);
	  dv[2]=0.5*(mx[nx*(j+0)+(i+1)]-mx[nx*(j+0)+(i-1)]);
	  dv[3]=0.5*(mx[nx*(j+1)+(i+1)]-mx[nx*(j+1)+(i-1)]);
	  dvxdx=func_bc(dv)*idx;

	  dv[0]=0.5*(my[nx*(j-2)+(i+1)]-my[nx*(j-2)+(i-1)]);
	  dv[1]=0.5*(my[nx*(j-1)+(i+1)]-my[nx*(j-1)+(i-1)]);
	  dv[2]=0.5*(my[nx*(j+0)+(i+1)]-my[nx*(j+0)+(i-1)]);
	  dv[3]=0.5*(my[nx*(j+1)+(i+1)]-my[nx*(j+1)+(i-1)]);
	  dvydx=func_bc(dv)*idx;

	  fy[4*ss+0]=-rnu*(dvxdy+dvydx);
	  fy[4*ss+1]=-rnu*(c4*dvydy-c2*dvxdx);
	  fy[4*ss+2]=-rnu*(dvzdy);
	  fy[4*ss+3]=fy[4*ss+0]*vx0+fy[4*ss+1]*vy0+fy[4*ss+2]*vz0;
	}
      }

      /* Update */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=yoff;j<ny-yoff;j++){
	int sip,sjp,si0,sj0;
	int sim,sip2,sjm,sjp2;
#pragma simd
	for (i=xoff;i<nx-xoff;i++){
	  ss=nx*j+i;
	  sip=nx*j+(i+1);
	  sjp=nx*(j+1)+i;
	  si0=nx*j+i;
	  sj0=nx*j+i;
	  sim=nx*j+(i-1);
	  sjm=nx*(j-1)+i;
	  sip2=nx*j+(i+2);
	  sjp2=nx*(j+2)+i;
	  double du[2][4];
	  /* dF/dx */
	  double valx[4][4]={{fx[4*sim+0],fx[4*si0+0],fx[4*sip+0],fx[4*sip2+0]},
			     {fx[4*sim+1],fx[4*si0+1],fx[4*sip+1],fx[4*sip2+1]},
			     {fx[4*sim+2],fx[4*si0+2],fx[4*sip+2],fx[4*sip2+2]},
			     {fx[4*sim+3],fx[4*si0+3],fx[4*sip+3],fx[4*sip2+3]}};
	  du[0][0]=func_df(&valx[0][0]);
	  du[0][1]=func_df(&valx[1][0]);
	  du[0][2]=func_df(&valx[2][0]);
	  du[0][3]=func_df(&valx[3][0]);
	  /* dG/dy */
	  double valy[4][4]={{fy[4*sjm+0],fy[4*sj0+0],fy[4*sjp+0],fy[4*sjp2+0]},
			     {fy[4*sjm+1],fy[4*sj0+1],fy[4*sjp+1],fy[4*sjp2+1]},
			     {fy[4*sjm+2],fy[4*sj0+2],fy[4*sjp+2],fy[4*sjp2+2]},
			     {fy[4*sjm+3],fy[4*sj0+3],fy[4*sjp+3],fy[4*sjp2+3]}};
	  du[1][0]=func_df(&valy[0][0]);
	  du[1][1]=func_df(&valy[1][0]);
	  du[1][2]=func_df(&valy[2][0]);
	  du[1][3]=func_df(&valy[3][0]);

	  mx[ss]=rk_fac[rk][0]*ut[4*ss+0]+rk_fac[rk][1]*(mx[ss]-(dtdx*du[0][0]+dtdy*du[1][0])*ir[ss]);
	  my[ss]=rk_fac[rk][0]*ut[4*ss+1]+rk_fac[rk][1]*(my[ss]-(dtdx*du[0][1]+dtdy*du[1][1])*ir[ss]);
	  mz[ss]=rk_fac[rk][0]*ut[4*ss+2]+rk_fac[rk][1]*(mz[ss]-(dtdx*du[0][2]+dtdy*du[1][2])*ir[ss]);
	  en[ss]=rk_fac[rk][0]*ut[4*ss+3]+rk_fac[rk][1]*(en[ss]-(dtdx*du[0][3]+dtdy*du[1][3]));
	}
      }

    } /* OpenMp */

    p[0]=mx;
    p[1]=my;
    p[2]=mz;
    mpi_sdrv2d(p,3,nx,ny,xoff,yoff,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(mx,nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(my,nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(mz,nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(mx,nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(my,nx,ny,xoff,yoff,0,-1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(mz,nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
  }
  p[0]=en;
  mpi_sdrv2d(p,1,nx,ny,xoff,yoff,mpi_rank,mpi_numx,mpi_numy);
  mpi_xbc2d(en,nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
  mpi_ybc2d(en,nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
    
  /* Return */
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      ss=nx*j+i;
      mx[ss]*=ro[ss];		/* vx => mx */
      my[ss]*=ro[ss];		/* vy => my */
      mz[ss]*=ro[ss];		/* vz => mz */
    }
  }
  
  free(ir);
  free(ut);
  free(fx);
  free(fy);
}

void hall_fd_ct_2d(double *ro, double *mx, double *my, double *mz,
		   double *en, double *bx, double *by, double *bz,
		   double dt, double dx, double dy,
		   int nx, int ny, int xoff, int yoff, double eta_h,
		   int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D finite-difference code for Hall term */
/* 2nd, 3rd, and 4th order */
/* LFF + Central-Upwind-CT */
/* Need four offsets (maximum) for boundary */
    
/* ro: density */
/* m?: moment */
/* en: energy normalized by B^2/(4 pi) */
/* b?: magnetic field */
/* eta_h: Hall resistivity (= va^2/(ion gyro freq.) */
{
  int i,j,ss,rk;
  const int ns=5;		/* number of grids in stencil */
  const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const int nxy=nx*ny;
  const double dtdx=dt/dx,dtdy=dt/dy;
  const double idx=1./dx,idy=1./dy;
  const double pi=4.0*atan(1.0);
  const double vphix=eta_h*pi*idx,vphiy=eta_h*pi*idy; /* Maximum whistler phase vel. */
  const double alpha=0.5;	/* Factor for hyper diffusion (<1) */
  double *ut,*ul,*ur,*ql,*qr,*fx,*fy;
  double *cx,*cy,*ez;
  double *ctx,*fc,eps=1e-6;
  double *hx,*hy,*hz;		/* Hall velocity, -j/nec */
  double *b2;
  double *d2u,*d2ul,*d2ur;			/* 2nd derivative for Hyper diffusion */
  double (*func_bc)(const double*)=&calc_bc;
  void (*func_lr)(const double*, double*, double*)=&calc_flr;
  void (*lfun_lr)(const double*, double*, double*)=&lcal_flr;
  double (*func_df)(const double*)=&calc_df;
  double (*func_dfc)(const double*)=&calc_dfc;
  double (*func_d2fc)(const double*)=&calc_d2fc;
  double *p[4];
  
  ut=(double*)malloc(sizeof(double)*4*nxy);
  ul=(double*)malloc(sizeof(double)*8*nxy);
  ur=(double*)malloc(sizeof(double)*8*nxy);
  ql=(double*)malloc(sizeof(double)*nxy);
  qr=(double*)malloc(sizeof(double)*nxy);
  fx=(double*)malloc(sizeof(double)*4*nxy);
  fy=(double*)malloc(sizeof(double)*4*nxy);
  cx=(double*)malloc(sizeof(double)*nxy);
  cy=(double*)malloc(sizeof(double)*nxy);
  ez=(double*)malloc(sizeof(double)*nxy);
  ctx=(double*)malloc(sizeof(double)*nxy);
  fc=(double*)malloc(sizeof(double)*nxy);
  hx=(double*)malloc(sizeof(double)*nxy);
  hy=(double*)malloc(sizeof(double)*nxy);
  hz=(double*)malloc(sizeof(double)*nxy);
  b2=(double*)malloc(sizeof(double)*nxy);
  d2u=(double*)malloc(sizeof(double)*3*nxy);
  d2ul=(double*)malloc(sizeof(double)*3*nxy);
  d2ur=(double*)malloc(sizeof(double)*3*nxy);

  /* Copy */
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      ss=nx*j+i;
      ut[4*ss+0]=bx[ss];
      ut[4*ss+1]=by[ss];
      ut[4*ss+2]=bz[ss];
      ut[4*ss+3]=en[ss];
    }
  }

  for (rk=0;rk<R_K;rk++){
    
#ifdef _OPENMP
#pragma omp parallel private(i,j,ss)
#endif
    {

      /* Calculate B at cell center */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=0;j<ny;j++){
#pragma simd
	for (i=1;i<nx-2;i++){
	  double val[]={bx[nx*j+(i-1)],bx[nx*j+(i+0)],
			bx[nx*j+(i+1)],bx[nx*j+(i+2)]};
	  cx[nx*j+i]=func_bc(val);
	}
      }
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=1;j<ny-2;j++){
#pragma simd
	for (i=0;i<nx;i++){
	  double val[]={by[nx*(j-1)+i],by[nx*(j+0)+i],
			by[nx*(j+1)+i],by[nx*(j+2)+i]};
	  cy[nx*j+i]=func_bc(val);
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	p[0]=cx;
	p[1]=cy;
	mpi_sdrv2d(p,2,nx,ny,xoff,yoff,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(&cx[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(&cy[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(&cx[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(&cy[0],nx,ny,xoff,yoff,0,-1,mpi_rank,mpi_numx,mpi_numy);
      }

      /* Calculate Hall velocity at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-2;j++){
	int sim2,sim1,sip0,sip1,sip2;
	int sjm2,sjm1,sjp0,sjp1,sjp2;
	double jx,jy,jz,tmp;
#pragma simd
	for (i=2;i<nx-2;i++){
	  ss=nx*j+i;
	  sim2=nx*j+(i-2);
	  sim1=nx*j+(i-1);
	  sip0=nx*j+(i+0);
	  sip1=nx*j+(i+1);
	  sip2=nx*j+(i+2);
	  sjm2=nx*(j-2)+i;
	  sjm1=nx*(j-1)+i;
	  sjp0=nx*(j+0)+i;
	  sjp1=nx*(j+1)+i;
	  sjp2=nx*(j+2)+i;
	  double dbx1[]={cx[sjm2],cx[sjm1],cx[sjp0],cx[sjp1],cx[sjp2]};
	  double dby0[]={cy[sim2],cy[sim1],cy[sip0],cy[sip1],cy[sip2]};
	  double dbz0[]={bz[sim2],bz[sim1],bz[sip0],bz[sip1],bz[sip2]};
	  double dbz1[]={bz[sjm2],bz[sjm1],bz[sjp0],bz[sjp1],bz[sjp2]};
	  jx=+func_dfc(dbz1)*idy;
	  jy=-func_dfc(dbz0)*idx;
	  jz=+func_dfc(dby0)*idx-func_dfc(dbx1)*idy;
	  tmp=eta_h/ro[ss];
	  hx[ss]=-tmp*jx;
	  hy[ss]=-tmp*jy;
	  hz[ss]=-tmp*jz;
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	p[0]=hx;
	p[1]=hy;
	p[2]=hz;
	mpi_sdrv2d(p,3,nx,ny,xoff,yoff,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(&hx[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(&hy[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(&hz[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(&hx[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(&hy[0],nx,ny,xoff,yoff,0,-1,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(&hz[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
      }

      /* Some variables */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
#pragma simd
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  ez[ss]=0.0;
	  ctx[ss]=0.5;
	  b2[ss]=0.5*(cx[ss]*cx[ss]+cy[ss]*cy[ss]+bz[ss]*bz[ss]);
	}
      }
      /* CT E-field 2D upwind weighting */
#if (CTW)
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=1;j<ny;j++){
	double vxc,vyc,roc,bxc,byc;
#pragma simd
	for (i=1;i<nx;i++){
	  /* Upwinding with respect to whistler mode */
	  vxc=0.25*(hx[nx*(j-1)+(i-1)]+hx[nx*(j-1)+i]+
		    hx[nx*j+(i-1)]+hx[nx*j+i]);
	  vyc=0.25*(hy[nx*(j-1)+(i-1)]+hy[nx*(j-1)+i]+
		    hy[nx*j+(i-1)]+hy[nx*j+i]);
	  roc=0.25*(ro[nx*(j-1)+(i-1)]+ro[nx*(j-1)+i]+
		    ro[nx*j+(i-1)]+ro[nx*j+i]);
	  roc=1.0/roc;
	  bxc=0.5*(bx[nx*(j-1)+i]+bx[nx*j+i]);
	  byc=0.5*(by[nx*j+(i-1)]+by[nx*j+i]);
	  vxc=fabs(vxc)+fabs(bxc*roc)*vphix;
	  vyc=fabs(vyc)+fabs(byc*roc)*vphiy;

	  ctx[nx*j+i]=(vxc+0.5*eps)/(vxc+vyc+eps);
	}
      }
#endif

      /* 2nd derivative in X for hyper diffusion */
#if (HDF)
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
	int sm2,sm,sp,sp2;
#pragma simd
	for (i=2;i<nx-2;i++){
	  ss=nx*j+i;
	  sm2=nx*j+(i-2);
	  sm=nx*j+(i-1);
	  sp=nx*j+(i+1);
	  sp2=nx*j+(i+2);
	  double val[3][5]={{cy[sm2],cy[sm],cy[ss],cy[sp],cy[sp2]},
	  		    {bz[sm2],bz[sm],bz[ss],bz[sp],bz[sp2]},
	  		    {b2[sm2],b2[sm],b2[ss],b2[sp],b2[sp2]}};

	  d2u[0*nxy+ss]=func_d2fc(&val[0][0]); /* by */
	  d2u[1*nxy+ss]=func_d2fc(&val[1][0]); /* bz */
	  d2u[2*nxy+ss]=func_d2fc(&val[2][0]); /* en */
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	p[0]=&d2u[0*nxy];
	p[1]=&d2u[1*nxy];
	p[2]=&d2u[2*nxy];
	mpi_sdrv2d(p,3,nx,ny,xoff,yoff,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(p[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(p[1],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(p[2],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(p[0],nx,ny,xoff,yoff,0,-1,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(p[1],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(p[2],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
      }
#endif

      /* Variables at cell face along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
	int sl,sr;
	int sm2,sm,sp,sp2;
#pragma simd
	for (i=2;i<nx-2;i++){
	  ss=nx*j+i;
	  sl=nx*j+(i+1);
	  sr=nx*j+i;
	  sm2=nx*j+(i-2);
	  sm=nx*j+(i-1);
	  sp=nx*j+(i+1);
	  sp2=nx*j+(i+2);
	  /* Variables in the stencil */
	  double ros[]={ro[sm2],ro[sm],ro[ss],ro[sp],ro[sp2]};
	  double hxs[]={hx[sm2],hx[sm],hx[ss],hx[sp],hx[sp2]};
	  double hys[]={hy[sm2],hy[sm],hy[ss],hy[sp],hy[sp2]};
	  double hzs[]={hz[sm2],hz[sm],hz[ss],hz[sp],hz[sp2]};
	  double bxs[]={cx[sm2],cx[sm],cx[ss],cx[sp],cx[sp2]};
	  double bys[]={cy[sm2],cy[sm],cy[ss],cy[sp],cy[sp2]};
	  double bzs[]={bz[sm2],bz[sm],bz[ss],bz[sp],bz[sp2]};
	  double ens[]={en[sm2],en[sm],en[ss],en[sp],en[sp2]};

	  double vl[7],vr[7];
	  func_lr(ros,&vl[0],&vr[0]);
	  func_lr(hxs,&vl[1],&vr[1]);
	  func_lr(hys,&vl[2],&vr[2]);
	  func_lr(hzs,&vl[3],&vr[3]);
	  func_lr(bys,&vl[4],&vr[4]);
	  func_lr(bzs,&vl[5],&vr[5]);
	  func_lr(ens,&vl[6],&vr[6]);
	  
	  /* Left-face @ i+1/2 */
	  ul[8*sl+0]=vl[0];	/* ro */
	  ul[8*sl+1]=vl[1];	/* hx */
	  ul[8*sl+2]=vl[2];	/* hy */
	  ul[8*sl+3]=vl[3];	/* hz */
	  ul[8*sl+4]=bx[sl];	/* bx */
	  ul[8*sl+5]=vl[4];	/* by */
	  ul[8*sl+6]=vl[5];	/* bz */
	  ul[8*sl+7]=vl[6];	/* en */
	  /* Right-face @ i-1/2 */
	  ur[8*sr+0]=vr[0];	/* ro */
	  ur[8*sr+1]=vr[1];	/* hx */
	  ur[8*sr+2]=vr[2];	/* hy */
	  ur[8*sr+3]=vr[3];	/* hz */
	  ur[8*sr+4]=bx[sr];	/* bx */
	  ur[8*sr+5]=vr[4];	/* by */
	  ur[8*sr+6]=vr[5];	/* bz */
	  ur[8*sr+7]=vr[6];	/* en */

	  /* Linear interpolation of numerical flux of By */
	  double val[]={cy[sm2]*hx[sm2],cy[sm]*hx[sm],cy[ss]*hx[ss],cy[sp]*hx[sp],cy[sp2]*hx[sp2]};
	  lfun_lr(val,&vl[0],&vr[0]);
	  lfun_lr(hys,&vl[1],&vr[1]);
	  ql[sl]=vl[0]-bx[sl]*vl[1];
	  qr[sr]=vr[0]-bx[sr]*vr[1];

#if (HDF)
	  /* Interpolation for hyper diffusion */
	  double d2us[3][5]={{d2u[0*nxy+sm2],d2u[0*nxy+sm],d2u[0*nxy+ss],d2u[0*nxy+sp],d2u[0*nxy+sp2]},
			     {d2u[1*nxy+sm2],d2u[1*nxy+sm],d2u[1*nxy+ss],d2u[1*nxy+sp],d2u[1*nxy+sp2]},
			     {d2u[2*nxy+sm2],d2u[2*nxy+sm],d2u[2*nxy+ss],d2u[2*nxy+sp],d2u[2*nxy+sp2]}};
	  func_lr(&d2us[0][0],&vl[0],&vr[0]);
	  func_lr(&d2us[1][0],&vl[1],&vr[1]);
	  func_lr(&d2us[2][0],&vl[2],&vr[2]);
	  d2ul[3*sl+0]=vl[0];	/* by */
	  d2ul[3*sl+1]=vl[1];	/* bz */
	  d2ul[3*sl+2]=vl[2];	/* en */
	  d2ur[3*sr+0]=vr[0];	/* by */
	  d2ur[3*sr+1]=vr[1];	/* bz */
	  d2ur[3*sr+2]=vr[2];	/* en */
#endif
	}
      }
      /* Numerical flux at cell face along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
#pragma simd
	for (i=3;i<nx-2;i++){
	  double bn,smax=0,flux[4];
	  ss=nx*j+i;
	  bn=0.5*(ul[8*ss+4]+ur[8*ss+4]);
	  smax+=max(fabs(ul[8*ss+1]),fabs(ur[8*ss+1])); /* Hall velocity */
	  smax+=vphix*fabs(bn/(0.5*(ul[8*ss+0]+ur[8*ss+0]))); /* Whistler velocity */
	  hall_flux_lf(ul[8*ss+0],ul[8*ss+1],ul[8*ss+2],ul[8*ss+3],ul[8*ss+5],ul[8*ss+6],ul[8*ss+7],
		       ur[8*ss+0],ur[8*ss+1],ur[8*ss+2],ur[8*ss+3],ur[8*ss+5],ur[8*ss+6],ur[8*ss+7],
		       bn,smax,
		       &flux[1],&flux[2],&flux[3]);
	  fx[4*ss+0]=0;		/* bx */
	  fx[4*ss+1]=flux[1];	/* by */
	  fx[4*ss+2]=flux[2];	/* bz */
	  fx[4*ss+3]=flux[3];	/* en */

#if (HDF)
	  /* Hyper diffusion term */
	  fx[4*ss+1]+=0.125*alpha*smax*(d2ur[3*ss+0]-d2ul[3*ss+0]);
	  fx[4*ss+2]+=0.125*alpha*smax*(d2ur[3*ss+1]-d2ul[3*ss+1]);
	  fx[4*ss+3]+=0.125*alpha*smax*(d2ur[3*ss+2]-d2ul[3*ss+2]);
#endif
	  
	  /* Divide central and upwind parts in numerical flux of By */
	  fc[ss]=0.5*(ql[ss]+qr[ss]); /* Central part */
	  fx[4*ss+1]-=fc[ss];	      /* Upwind part */
	}
      }
      /* Numerical flux of By at cell corner along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-2;j++){
	int sl,sr;
	int sm2,sm,sp,sp2;
#pragma simd
	for (i=3;i<nx-2;i++){
	  ss=nx*j+i;
	  sl=nx*(j+1)+i;
	  sr=nx*j+i;
	  sm2=nx*(j-2)+i;
	  sm=nx*(j-1)+i;
	  sp=nx*(j+1)+i;
	  sp2=nx*(j+2)+i;
	  double val[]={fx[4*sm2+1],fx[4*sm+1],fx[4*ss+1],fx[4*sp+1],fx[4*sp2+1]};
	  double vac[]={fc[sm2],fc[sm],fc[ss],fc[sp],fc[sp2]};

	  double vl[2],vr[2];
	  lfun_lr(val,&vl[0],&vr[0]); /* Upwind part */
	  lfun_lr(vac,&vl[1],&vr[1]); /* Central part */
	  ul[2*sl+0]=vl[0];
	  ur[2*sr+0]=vr[0];
	  ul[2*sl+1]=vl[1];
	  ur[2*sr+1]=vr[1];
	}
      }
      /* E-field at cell corner along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
#pragma simd
	for (i=3;i<nx-2;i++){
	  ss=nx*j+i;
	  ez[ss]+=-0.5*((ul[2*ss+0]+ur[2*ss+0])+(1.0-ctx[ss])*(ul[2*ss+1]+ur[2*ss+1]));
	}
      }

      /* 2nd derivative in Y for hyper diffusion */
#if (HDF)
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-2;j++){
	int sm2,sm,sp,sp2;
#pragma simd
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  sm2=nx*(j-2)+i;
	  sm=nx*(j-1)+i;
	  sp=nx*(j+1)+i;
	  sp2=nx*(j+2)+i;
	  double val[3][5]={{bz[sm2],bz[sm],bz[ss],bz[sp],bz[sp2]},
	  		    {cx[sm2],cx[sm],cx[ss],cx[sp],cx[sp2]},
	  		    {b2[sm2],b2[sm],b2[ss],b2[sp],b2[sp2]}};
	  
	  d2u[0*nxy+ss]=func_d2fc(&val[0][0]); /* bz */
	  d2u[1*nxy+ss]=func_d2fc(&val[1][0]); /* bx */
	  d2u[2*nxy+ss]=func_d2fc(&val[2][0]); /* en */
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	p[0]=&d2u[0*nxy];
	p[1]=&d2u[1*nxy];
	p[2]=&d2u[2*nxy];
	mpi_sdrv2d(p,3,nx,ny,xoff,yoff,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(p[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(p[1],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_xbc2d(p[2],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(p[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(p[1],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(p[2],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
      }
#endif
      
      /* Variables at cell face along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-2;j++){
	int sl,sr;
	int sm2,sm,sp,sp2;
#pragma simd
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  sl=nx*(j+1)+i;
	  sr=nx*j+i;
	  sm2=nx*(j-2)+i;
	  sm=nx*(j-1)+i;
	  sp=nx*(j+1)+i;
	  sp2=nx*(j+2)+i;
	  /* Variables in the stencil */
	  double ros[]={ro[sm2],ro[sm],ro[ss],ro[sp],ro[sp2]};
	  double hxs[]={hx[sm2],hx[sm],hx[ss],hx[sp],hx[sp2]};
	  double hys[]={hy[sm2],hy[sm],hy[ss],hy[sp],hy[sp2]};
	  double hzs[]={hz[sm2],hz[sm],hz[ss],hz[sp],hz[sp2]};
	  double bxs[]={cx[sm2],cx[sm],cx[ss],cx[sp],cx[sp2]};
	  double bys[]={cy[sm2],cy[sm],cy[ss],cy[sp],cy[sp2]};
	  double bzs[]={bz[sm2],bz[sm],bz[ss],bz[sp],bz[sp2]};
	  double ens[]={en[sm2],en[sm],en[ss],en[sp],en[sp2]};

	  double vl[7],vr[7];
	  func_lr(ros,&vl[0],&vr[0]);
	  func_lr(hys,&vl[1],&vr[1]);
	  func_lr(hzs,&vl[2],&vr[2]);
	  func_lr(hxs,&vl[3],&vr[3]);	  
	  func_lr(bzs,&vl[4],&vr[4]);
	  func_lr(bxs,&vl[5],&vr[5]);	  
	  func_lr(ens,&vl[6],&vr[6]);

	  /* Left-face @ j+1/2 */
	  ul[8*sl+0]=vl[0];	/* ro */
	  ul[8*sl+1]=vl[1];	/* hy */
	  ul[8*sl+2]=vl[2];	/* hz */
	  ul[8*sl+3]=vl[3];	/* hx */
	  ul[8*sl+4]=by[sl];	/* by */
	  ul[8*sl+5]=vl[4];	/* bz */
	  ul[8*sl+6]=vl[5];	/* bx */
	  ul[8*sl+7]=vl[6];	/* en */
	  /* Right-face @ j-1/2 */
	  ur[8*sr+0]=vr[0];	/* ro */
	  ur[8*sr+1]=vr[1];	/* hy */
	  ur[8*sr+2]=vr[2];	/* hz */
	  ur[8*sr+3]=vr[3];	/* hx */
	  ur[8*sr+4]=by[sr];	/* by */
	  ur[8*sr+5]=vr[4];	/* bz */
	  ur[8*sr+6]=vr[5];	/* bx */
	  ur[8*sr+7]=vr[6];	/* en */

	  /* Linear interpolation of numerical flux of Bx */
	  double val[]={cx[sm2]*hy[sm2],cx[sm]*hy[sm],cx[ss]*hy[ss],cx[sp]*hy[sp],cx[sp2]*hy[sp2]};
	  lfun_lr(val,&vl[0],&vr[0]);
	  lfun_lr(hxs,&vl[1],&vr[1]);
	  ql[sl]=vl[0]-by[sl]*vl[1];
	  qr[sr]=vr[0]-by[sr]*vr[1];

#if (HDF)
	  /* Interpolation for hyper diffusion */
	  double d2us[3][5]={{d2u[0*nxy+sm2],d2u[0*nxy+sm],d2u[0*nxy+ss],d2u[0*nxy+sp],d2u[0*nxy+sp2]},
			     {d2u[1*nxy+sm2],d2u[1*nxy+sm],d2u[1*nxy+ss],d2u[1*nxy+sp],d2u[1*nxy+sp2]},
			     {d2u[2*nxy+sm2],d2u[2*nxy+sm],d2u[2*nxy+ss],d2u[2*nxy+sp],d2u[2*nxy+sp2]}};
	  func_lr(&d2us[0][0],&vl[0],&vr[0]);
	  func_lr(&d2us[1][0],&vl[1],&vr[1]);
	  func_lr(&d2us[2][0],&vl[2],&vr[2]);
	  d2ul[3*sl+0]=vl[0];	/* bz */
	  d2ul[3*sl+1]=vl[1];	/* bx */
	  d2ul[3*sl+2]=vl[2];	/* en */
	  d2ur[3*sr+0]=vr[0];	/* bz */
	  d2ur[3*sr+1]=vr[1];	/* bx */
	  d2ur[3*sr+2]=vr[2];	/* en */
#endif
	}
      }
      /* Numerical flux at cell face along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
#pragma simd
	for (i=0;i<nx;i++){
	  double bn,smax=0,flux[4];
	  ss=nx*j+i;
	  bn=0.5*(ul[8*ss+4]+ur[8*ss+4]);
	  smax+=max(fabs(ul[8*ss+1]),fabs(ur[8*ss+1])); /* Hall velocity */
	  smax+=vphiy*fabs(bn/(0.5*(ul[8*ss+0]+ur[8*ss+0]))); /* Whistler velocity */
	  hall_flux_lf(ul[8*ss+0],ul[8*ss+1],ul[8*ss+2],ul[8*ss+3],ul[8*ss+5],ul[8*ss+6],ul[8*ss+7],
		       ur[8*ss+0],ur[8*ss+1],ur[8*ss+2],ur[8*ss+3],ur[8*ss+5],ur[8*ss+6],ur[8*ss+7],
		       bn,smax,
		       &flux[1],&flux[2],&flux[3]);
	  fy[4*ss+1]=0;		/* by */
	  fy[4*ss+2]=flux[1];	/* bz */
	  fy[4*ss+0]=flux[2];	/* bx */
	  fy[4*ss+3]=flux[3];	/* en */

#if (HDF)
	  /* Hyper diffusion term */
	  fy[4*ss+2]+=0.125*alpha*smax*(d2ur[3*ss+0]-d2ul[3*ss+0]);
	  fy[4*ss+0]+=0.125*alpha*smax*(d2ur[3*ss+1]-d2ul[3*ss+1]);
	  fy[4*ss+3]+=0.125*alpha*smax*(d2ur[3*ss+2]-d2ul[3*ss+2]);
#endif

	  /* Divide central and upwind parts in numerical flux of Bx */
	  fc[ss]=0.5*(ql[ss]+qr[ss]); /* Central part */
	  fy[4*ss+0]-=fc[ss];	      /* Upwind part */
	}
      }
      /* Numerical flux of Bx at cell corner along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
	int sl,sr;
	int sm2,sm,sp,sp2;
#pragma simd
	for (i=2;i<nx-2;i++){
	  ss=nx*j+i;
	  sl=nx*j+(i+1);
	  sr=nx*j+i;
	  sm2=nx*j+(i-2);
	  sm=nx*j+(i-1);
	  sp=nx*j+(i+1);
	  sp2=nx*j+(i+2);
	  double val[]={fy[4*sm2+0],fy[4*sm+0],fy[4*ss+0],fy[4*sp+0],fy[4*sp2+0]};
	  double vac[]={fc[sm2],fc[sm],fc[ss],fc[sp],fc[sp2]};

	  double vl[2],vr[2];
	  lfun_lr(val,&vl[0],&vr[0]); /* Upwind part */
	  lfun_lr(vac,&vl[1],&vr[1]); /* Central part */
	  ul[2*sl+0]=vl[0];
	  ur[2*sr+0]=vr[0];
	  ul[2*sl+1]=vl[1];
	  ur[2*sr+1]=vr[1];
	}
      }
      /* E-field at cell corner along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
#pragma simd
	for (i=3;i<nx-2;i++){
	  ss=nx*j+i;
	  ez[ss]+=+0.5*((ul[2*ss+0]+ur[2*ss+0])+ctx[ss]*(ul[2*ss+1]+ur[2*ss+1]));
	}
      }

      /* Update variable at cell center */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=yoff;j<ny-yoff;j++){
	int sip,sjp,si0,sj0;
	int sim,sip2,sjm,sjp2;
#pragma simd
	for (i=xoff;i<nx-xoff;i++){
	  ss=nx*j+i;
	  sip=nx*j+(i+1);
	  sjp=nx*(j+1)+i;
	  si0=nx*j+i;
	  sj0=nx*j+i;
	  sim=nx*j+(i-1);
	  sjm=nx*(j-1)+i;
	  sip2=nx*j+(i+2);
	  sjp2=nx*(j+2)+i;
	  double du[2][4];
	  /* dF/dx */
	  double valx[4][4]={{fx[4*sim+0],fx[4*si0+0],fx[4*sip+0],fx[4*sip2+0]},
			     {fx[4*sim+1],fx[4*si0+1],fx[4*sip+1],fx[4*sip2+1]},
			     {fx[4*sim+2],fx[4*si0+2],fx[4*sip+2],fx[4*sip2+2]},
			     {fx[4*sim+3],fx[4*si0+3],fx[4*sip+3],fx[4*sip2+3]}};
	  du[0][2]=func_df(&valx[2][0]);
	  du[0][3]=func_df(&valx[3][0]);
	  /* dG/dy */
	  double valy[4][4]={{fy[4*sjm+0],fy[4*sj0+0],fy[4*sjp+0],fy[4*sjp2+0]},
			     {fy[4*sjm+1],fy[4*sj0+1],fy[4*sjp+1],fy[4*sjp2+1]},
			     {fy[4*sjm+2],fy[4*sj0+2],fy[4*sjp+2],fy[4*sjp2+2]},
			     {fy[4*sjm+3],fy[4*sj0+3],fy[4*sjp+3],fy[4*sjp2+3]}};
	  du[1][2]=func_df(&valy[2][0]);
	  du[1][3]=func_df(&valy[3][0]);

	  bz[ss]=rk_fac[rk][0]*ut[4*ss+2]+rk_fac[rk][1]*(bz[ss]-dtdx*du[0][2]-dtdy*du[1][2]);
	  en[ss]=rk_fac[rk][0]*ut[4*ss+3]+rk_fac[rk][1]*(en[ss]-dtdx*du[0][3]-dtdy*du[1][3]);
	}
      }
      /* Update CT Bx */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=yoff;j<ny-yoff;j++){
	int sp,s0;
	int sm,sp2;
	double du;
#pragma simd
	for (i=xoff;i<nx-xoff+1;i++){
	  ss=nx*j+i;
	  sp=nx*(j+1)+i;
	  s0=nx*j+i;
	  sm=nx*(j-1)+i;
	  sp2=nx*(j+2)+i;
	  double val[]={ez[sm],ez[s0],ez[sp],ez[sp2]};
	  du=func_df(val);
	  bx[ss]=rk_fac[rk][0]*ut[4*ss+0]+rk_fac[rk][1]*(bx[ss]-dtdy*du);
	}
      }
      /* Update CT By */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=yoff;j<ny-yoff+1;j++){
	int sp,s0;
	int sm,sp2;
	double du;
#pragma simd
	for (i=xoff;i<nx-xoff;i++){
	  ss=nx*j+i;
	  sp=nx*j+(i+1);
	  s0=nx*j+i;
	  sm=nx*j+(i-1);
	  sp2=nx*j+(i+2);
	  double val[]={ez[sm],ez[s0],ez[sp],ez[sp2]};
	  du=func_df(val);
	  by[ss]=rk_fac[rk][0]*ut[4*ss+1]+rk_fac[rk][1]*(by[ss]+dtdx*du);
	}
      }

    } /* OpenMP */
    p[0]=en;
    p[1]=bx;
    p[2]=by;
    p[3]=bz;
    mpi_sdrv2d(p,4,nx,ny,xoff,yoff,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(&en[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(&bx[0],nx,ny,xoff,yoff,1,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(&by[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_xbc2d(&bz[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&en[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&bx[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&by[0],nx,ny,xoff,yoff,1,-1,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&bz[0],nx,ny,xoff,yoff,0,+1,mpi_rank,mpi_numx,mpi_numy);
  }

  free(ut);
  free(ul);
  free(ur);
  free(ql);
  free(qr);
  free(fx);
  free(fy);
  free(cx);
  free(cy);
  free(ez);
  free(ctx);
  free(fc);
  free(hx);
  free(hy);
  free(hz);
  free(b2);
  free(d2u);
  free(d2ul);
  free(d2ur);
}
