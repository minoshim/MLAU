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
  (*fl)=0.25*(-f[1]+4.0*f[2]+f[3]);
  (*fr)=0.25*(-f[3]+4.0*f[2]+f[1]);
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
  muscl_mm_cal_flr(f[1],f[2],f[3],fl,fr);
  /* muscl_mc_cal_flr(f[1],f[2],f[3],fl,fr); */
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
	mpi_ybc2d(&cx[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(&cy[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
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
	mpi_ybc2d(&dvx[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
	mpi_ybc2d(&dvy[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
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
    mpi_ybc2d(&ro[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&mx[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&my[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&mz[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&en[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&bx[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&by[0],nx,ny,xoff,yoff,1,+0,mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(&bz[0],nx,ny,xoff,yoff,0,+0,mpi_rank,mpi_numx,mpi_numy);
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
