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

void mpi_sdrv3d(double *f[], int nn, int nx, int ny, int nz, int xoff, int yoff, int zoff,
		int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);
void mpi_xbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);
void mpi_ybc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);
void mpi_zbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);

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

void mhd_fd_ct_3d(double *ro, double *mx, double *my, double *mz,
		  double *en, double *bx, double *by, double *bz,
		  double dt, double dx, double dy, double dz,
		  int nx, int ny, int nz, int xoff, int yoff, int zoff, double gamma,
		  int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
		  
/* 3D finite-difference code for MHD */
/* 2nd, 3rd, and 4th order with characteristic decomposition */
/* MLAU/HLLD/ROE + Central-Upwind-CT */
/* Need four offsets (maximum) for boundary */

/* ro: density */
/* m?: moment */
/* en: energy normalized by B^2/(4 pi) */
/* b?: magnetic field */
/* gamma: specific heat ratio */

{
  int i,j,k,ss,s1;
  int rk;
  const int ns=5;		/* number of grids in stencil */
  const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const int nxyz=nx*ny*nz;
  const double dtdx=dt/dx,dtdy=dt/dy,dtdz=dt/dz;
  double *ut,*fx,*fy,*fz;
  double *vx,*vy,*vz,*pr,*cx,*cy,*cz,*ex,*ey,*ez;
  double *ct,*fc,eps=1e-6;
  double *m2,*dvx,*dvy,*dvz;
  double (*func_bc)(const double*)=&calc_bc;
  void (*func_lr)(const double*, double*, double*)=&calc_flr;
  void (*lfun_lr)(const double*, double*, double*)=&lcal_flr;
  double (*func_df)(const double*)=&calc_df;
  void (*func_flux)(double, double, double, double, double, double, double,
		    double, double, double, double, double, double, double,
		    double, double, double*,
		    double*, double*, double*, double*, double*, double*, double*);
  /* func_flux=&calc_flux_mlau; */
  /* func_flux=&calc_flux_hlld; */
  /* func_flux=&calc_flux_roe; */
  func_flux=&calc_flux_lhlld;
  double *p[8];
  
  ut=(double*)malloc(sizeof(double)*8*nxyz);
  fx=(double*)malloc(sizeof(double)*8*nxyz);
  fy=(double*)malloc(sizeof(double)*8*nxyz);
  fz=(double*)malloc(sizeof(double)*8*nxyz);
  vx=(double*)malloc(sizeof(double)*nxyz);
  vy=(double*)malloc(sizeof(double)*nxyz);
  vz=(double*)malloc(sizeof(double)*nxyz);
  pr=(double*)malloc(sizeof(double)*nxyz);
  cx=(double*)malloc(sizeof(double)*nxyz);
  cy=(double*)malloc(sizeof(double)*nxyz);
  cz=(double*)malloc(sizeof(double)*nxyz);
  ex=(double*)malloc(sizeof(double)*nxyz);
  ey=(double*)malloc(sizeof(double)*nxyz);
  ez=(double*)malloc(sizeof(double)*nxyz);
  ct=(double*)malloc(sizeof(double)*3*nxyz);
  fc=(double*)malloc(sizeof(double)*2*nxyz);
  m2=(double*)malloc(sizeof(double)*nxyz);
  dvx=(double*)malloc(sizeof(double)*nxyz);
  dvy=(double*)malloc(sizeof(double)*nxyz);
  dvz=(double*)malloc(sizeof(double)*nxyz);

  /* Copy */
  for (k=0;k<nz;k++){
    for (j=0;j<ny;j++){
      for (i=0;i<nx;i++){
	ss=nx*(ny*k+j)+i;
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
  }
  
  for (rk=0;rk<R_K;rk++){
    
#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ss,s1)
#endif
    {

      /* Calculate B at cell center */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (k=0;k<nz;k++){
	for (j=0;j<ny;j++){
#pragma simd
	  for (i=1;i<nx-2;i++){
	    double val[]={bx[nx*(ny*k+j)+(i-1)],bx[nx*(ny*k+j)+(i+0)],
			  bx[nx*(ny*k+j)+(i+1)],bx[nx*(ny*k+j)+(i+2)]};
	    cx[nx*(ny*k+j)+i]=func_bc(val);
	  }
	}
	for (j=1;j<ny-2;j++){
#pragma simd
	  for (i=0;i<nx;i++){
	    double val[]={by[nx*(ny*k+(j-1))+i],by[nx*(ny*k+(j+0))+i],
			  by[nx*(ny*k+(j+1))+i],by[nx*(ny*k+(j+2))+i]};
	    cy[nx*(ny*k+j)+i]=func_bc(val);
	  }
	}
      }
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=1;k<nz-2;k++){
	for (j=0;j<ny;j++){
#pragma simd
	  for (i=0;i<nx;i++){
	    double val[]={bz[nx*(ny*(k-1)+j)+i],bz[nx*(ny*(k+0)+j)+i],
			  bz[nx*(ny*(k+1)+j)+i],bz[nx*(ny*(k+2)+j)+i]};
	    cz[nx*(ny*k+j)+i]=func_bc(val);
	  }
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	p[0]=cx;
	p[1]=cy;
	p[2]=cz;
	mpi_sdrv3d(p,3,nx,ny,nz,xoff,yoff,zoff,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_xbc3d(&cx[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_xbc3d(&cy[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_xbc3d(&cz[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_ybc3d(&cx[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_ybc3d(&cy[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_ybc3d(&cz[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_zbc3d(&cx[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_zbc3d(&cy[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_zbc3d(&cz[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
      }

      /* Primitive variable at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	double iro,v2,b2;
	for (j=0;j<ny;j++){
#pragma simd
	  for (i=0;i<nx;i++){
	    ss=nx*(ny*k+j)+i;
	    iro=1.0/ro[ss];
	    vx[ss]=mx[ss]*iro;
	    vy[ss]=my[ss]*iro;
	    vz[ss]=mz[ss]*iro;
	    v2=(vx[ss]*vx[ss]+vy[ss]*vy[ss]+vz[ss]*vz[ss]);
	    b2=(cx[ss]*cx[ss]+cy[ss]*cy[ss]+cz[ss]*cz[ss]);
	    pr[ss]=(gamma-1)*(en[ss]-0.5*(ro[ss]*v2+b2));
	    m2[ss]=v2/(gamma*pr[ss]*iro);
	    ex[ss]=ey[ss]=ez[ss]=0.0; /* Necessary initialize at cell corner */
	    ct[nxyz*0+ss]=ct[nxyz*1+ss]=ct[nxyz*2+ss]=0.5;
	  }
	}
      }

      /* CT E-field multi-D upwind weighting */
#if (CTW)
      /* X-Y plane for Ez */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (k=0;k<nz;k++){
	double vxc,vyc,roc,bxc,byc;
	for (j=1;j<ny;j++){
#pragma simd
	  for (i=1;i<nx;i++){
	    ss=nx*(ny*k+j)+i;
	    /* Upwinding with respect to Alfven mode */
	    vxc=0.25*(vx[nx*(ny*k+(j-1))+(i-1)]+vx[nx*(ny*k+(j-1))+i]+
		      vx[nx*(ny*k+j)+(i-1)]+vx[ss]);
	    vyc=0.25*(vy[nx*(ny*k+(j-1))+(i-1)]+vy[nx*(ny*k+(j-1))+i]+
		      vy[nx*(ny*k+j)+(i-1)]+vy[ss]);
	    roc=0.25*(ro[nx*(ny*k+(j-1))+(i-1)]+ro[nx*(ny*k+(j-1))+i]+
		      ro[nx*(ny*k+j)+(i-1)]+ro[ss]);
	    roc=1.0/sqrt(roc);
	    bxc=0.5*(bx[nx*(ny*k+(j-1))+i]+bx[ss]);
	    byc=0.5*(by[nx*(ny*k+j)+(i-1)]+by[ss]);
	    vxc=fabs(vxc)+fabs(bxc*roc);
	    vyc=fabs(vyc)+fabs(byc*roc);
	    
	    ct[nxyz*0+ss]=(vxc+0.5*eps)/(vxc+vyc+eps);
	  }
	}
      }
      /* Y-Z plane for Ex, Z-X plane for Ey */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=1;k<nz;k++){
	double vxc,vyc,vzc,roc,bxc,byc,bzc;
	for (j=1;j<ny;j++){
#pragma simd
	  for (i=0;i<nx;i++){
	    ss=nx*(ny*k+j)+i;
	    /* Upwinding with respect to Alfven mode */
	    vyc=0.25*(vy[nx*(ny*(k-1)+(j-1))+i]+vy[nx*(ny*(k-1)+j)+i]+
		      vy[nx*(ny*k+(j-1))+i]+vy[ss]);
	    vzc=0.25*(vz[nx*(ny*(k-1)+(j-1))+i]+vz[nx*(ny*(k-1)+j)+i]+
		      vz[nx*(ny*k+(j-1))+i]+vz[ss]);
	    roc=0.25*(ro[nx*(ny*(k-1)+(j-1))+i]+ro[nx*(ny*(k-1)+j)+i]+
		      ro[nx*(ny*k+(j-1))+i]+ro[ss]);
	    roc=1.0/sqrt(roc);
	    byc=0.5*(by[nx*(ny*(k-1)+j)+i]+by[ss]);
	    bzc=0.5*(bz[nx*(ny*k+(j-1))+i]+bz[ss]);
	    vyc=fabs(vyc)+fabs(byc*roc);
	    vzc=fabs(vzc)+fabs(bzc*roc);
	    
	    ct[nxyz*1+ss]=(vyc+0.5*eps)/(vyc+vzc+eps);
	  }
	}
	for (j=0;j<ny;j++){
#pragma simd
	  for (i=1;i<nx;i++){
	    ss=nx*(ny*k+j)+i;
	    /* Upwinding with respect to Alfven mode */
	    vzc=0.25*(vz[nx*(ny*(k-1)+j)+(i-1)]+vz[nx*(ny*(k-1)+j)+i]+
		      vz[nx*(ny*k+j)+(i-1)]+vz[ss]);
	    vxc=0.25*(vx[nx*(ny*(k-1)+j)+(i-1)]+vx[nx*(ny*(k-1)+j)+i]+
		      vx[nx*(ny*k+j)+(i-1)]+vx[ss]);
	    roc=0.25*(ro[nx*(ny*(k-1)+j)+(i-1)]+ro[nx*(ny*(k-1)+j)+i]+
		      ro[nx*(ny*k+j)+(i-1)]+ro[ss]);
	    roc=1.0/sqrt(roc);
	    bzc=0.5*(bz[nx*(ny*k+j)+(i-1)]+bz[ss]);
	    bxc=0.5*(bx[nx*(ny*(k-1)+j)+i]+bx[ss]);
	    vzc=fabs(vzc)+fabs(bzc*roc);
	    vxc=fabs(vxc)+fabs(bxc*roc);
	    
	    ct[nxyz*2+ss]=(vzc+0.5*eps)/(vzc+vxc+eps);
	  }
	}
      }
#endif	/* CTW */

      /* Calculate dvx,y,z for shock detection */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=1;k<nz-1;k++){
	for (j=1;j<ny-1;j++){
#pragma simd
	  for (i=1;i<nx-1;i++){
	    ss=nx*(ny*k+j)+i;
	    dvx[ss]=min((vx[nx*(ny*k+j)+(i+1)]-vx[nx*(ny*k+j)+i]),
			(vx[nx*(ny*k+j)+i]-vx[nx*(ny*k+j)+(i-1)]));
	    dvy[ss]=min((vy[nx*(ny*k+(j+1))+i]-vy[nx*(ny*k+j)+i]),
			(vy[nx*(ny*k+j)+i]-vy[nx*(ny*k+(j-1))+i]));
	    dvz[ss]=min((vz[nx*(ny*(k+1)+j)+i]-vz[nx*(ny*k+j)+i]),
			(vz[nx*(ny*k+j)+i]-vz[nx*(ny*(k-1)+j)+i]));
	  }
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	p[0]=dvx;
	p[1]=dvy;
	p[2]=dvz;
	mpi_sdrv3d(p,3,nx,ny,nz,xoff,yoff,zoff,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_xbc3d(&dvx[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_xbc3d(&dvy[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_xbc3d(&dvz[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_ybc3d(&dvx[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_ybc3d(&dvy[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_ybc3d(&dvz[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_zbc3d(&dvx[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_zbc3d(&dvy[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
	mpi_zbc3d(&dvz[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
      }
      
      /* Calculate numerical flux at cell face along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	int s1l,s1r;
	int sl,sr;
	int sm2,sm,sp,sp2;
	double *ul,*ur,*ql,*qr;
	ul=(double*)malloc(sizeof(double)*8*nx);
	ur=(double*)malloc(sizeof(double)*8*nx);
	ql=(double*)malloc(sizeof(double)*2*nx);
	qr=(double*)malloc(sizeof(double)*2*nx);
	for (j=0;j<ny;j++){
	  /* Interpolation */
	  for (i=2;i<nx-2;i++){
	    s1l=i+1;
	    s1r=i;
	    ss=nx*(ny*k+j)+i;
	    sl=nx*(ny*k+j)+(i+1);
	    sr=nx*(ny*k+j)+i;
	    sm2=nx*(ny*k+j)+(i-2);
	    sm=nx*(ny*k+j)+(i-1);
	    sp=nx*(ny*k+j)+(i+1);
	    sp2=nx*(ny*k+j)+(i+2);
	    /* Primitive variables in the stencil */
	    double ros[]={ro[sm2],ro[sm],ro[ss],ro[sp],ro[sp2]};
	    double vxs[]={vx[sm2],vx[sm],vx[ss],vx[sp],vx[sp2]};
	    double vys[]={vy[sm2],vy[sm],vy[ss],vy[sp],vy[sp2]};
	    double vzs[]={vz[sm2],vz[sm],vz[ss],vz[sp],vz[sp2]};
	    double bxs[]={cx[sm2],cx[sm],cx[ss],cx[sp],cx[sp2]};
	    double bys[]={cy[sm2],cy[sm],cy[ss],cy[sp],cy[sp2]};
	    double bzs[]={cz[sm2],cz[sm],cz[ss],cz[sp],cz[sp2]};
	    double prs[]={pr[sm2],pr[sm],pr[ss],pr[sp],pr[sp2]};

	    double vl[7],vr[7];

	    if (max(m2[sm],max(m2[ss],m2[sp])) >= 1.0){
	      mhd_c_reconst(ros,vxs,vys,vzs,bys,bzs,prs,bxs[ns/2],gamma,ns,vl,vr,func_lr);
	    } else{
	      mhd_a_reconst(ros,vxs,vys,vzs,bys,bzs,prs,bxs[ns/2],gamma,ns,vl,vr,func_lr);
	    }

	    /* Left-face @ i+1/2 */
	    ul[8*s1l+0]=vl[0];	/* ro */
	    ul[8*s1l+1]=vl[1];	/* vx */
	    ul[8*s1l+2]=vl[2];	/* vy */
	    ul[8*s1l+3]=vl[3];	/* vz */
	    ul[8*s1l+4]=bx[sl];	/* bx */
	    ul[8*s1l+5]=vl[4];	/* by */
	    ul[8*s1l+6]=vl[5];	/* bz */
	    ul[8*s1l+7]=vl[6];	/* pr */
	    /* Right-face @ i-1/2 */
	    ur[8*s1r+0]=vr[0];	/* ro */
	    ur[8*s1r+1]=vr[1];	/* vx */
	    ur[8*s1r+2]=vr[2];	/* vy */
	    ur[8*s1r+3]=vr[3];	/* vz */
	    ur[8*s1r+4]=bx[sr];	/* bx */
	    ur[8*s1r+5]=vr[4];	/* by */
	    ur[8*s1r+6]=vr[5];	/* bz */
	    ur[8*s1r+7]=vr[6];	/* pr */

	    /* Linear interpolation of numerical flux of By */
	    double vay[]={cy[sm2]*vx[sm2],cy[sm]*vx[sm],cy[ss]*vx[ss],cy[sp]*vx[sp],cy[sp2]*vx[sp2]};
	    lfun_lr(vay,&vl[0],&vr[0]);
	    lfun_lr(vys,&vl[1],&vr[1]);
	    ql[2*s1l+0]=vl[0]-bx[sl]*vl[1];
	    qr[2*s1r+0]=vr[0]-bx[sr]*vr[1];

	    /* Linear interpolation of numerical flux of Bz */
	    double vaz[]={cz[sm2]*vx[sm2],cz[sm]*vx[sm],cz[ss]*vx[ss],cz[sp]*vx[sp],cz[sp2]*vx[sp2]};
	    lfun_lr(vaz,&vl[0],&vr[0]);
	    lfun_lr(vzs,&vl[1],&vr[1]);
	    ql[2*s1l+1]=vl[0]-bx[sl]*vl[1];
	    qr[2*s1r+1]=vr[0]-bx[sr]*vr[1];
	  }
	  /* Riemann solver */
#pragma simd
	  for (i=3;i<nx-2;i++){
	    double bn,flux[8];
	    s1=i;
	    ss=nx*(ny*k+j)+i;
	    bn=0.5*(ul[8*s1+4]+ur[8*s1+4]);

	    double dvtm=min(min(dvy[nx*(ny*k+j)+(i-1)],dvy[nx*(ny*k+j)+i]),
			    min(dvz[nx*(ny*k+j)+(i-1)],dvz[nx*(ny*k+j)+i]));
	    double dvsd[2]={(vx[nx*(ny*k+j)+i]-vx[nx*(ny*k+j)+(i-1)]),dvtm};
	    func_flux(ul[8*s1+0],ul[8*s1+1],ul[8*s1+2],ul[8*s1+3],ul[8*s1+5],ul[8*s1+6],ul[8*s1+7],
		      ur[8*s1+0],ur[8*s1+1],ur[8*s1+2],ur[8*s1+3],ur[8*s1+5],ur[8*s1+6],ur[8*s1+7],
		      bn,gamma,dvsd,
		      &flux[0],&flux[1],&flux[2],&flux[3],&flux[5],&flux[6],&flux[7]);

	    fx[8*ss+0]=flux[0];	/* ro */
	    fx[8*ss+1]=flux[1];	/* mx */
	    fx[8*ss+2]=flux[2];	/* my */
	    fx[8*ss+3]=flux[3];	/* mz */
	    fx[8*ss+4]=0;	/* bx */
	    fx[8*ss+5]=flux[5];	/* by */
	    fx[8*ss+6]=flux[6];	/* bz */
	    fx[8*ss+7]=flux[7];	/* en */

	    /* Divide central and upwind parts in numerical flux of By */
	    fc[2*ss+0]=0.5*(ql[2*s1+0]+qr[2*s1+0]); /* Central part */
	    fx[8*ss+5]-=fc[2*ss+0];			     /* Upwind part */

	    /* Divide central and upwind parts in numerical flux of Bz */
	    fc[2*ss+1]=0.5*(ql[2*s1+1]+qr[2*s1+1]); /* Central part */
	    fx[8*ss+6]-=fc[2*ss+1];			     /* Upwind part */
	  }
	}
	free(ul);
	free(ur);
	free(ql);
	free(qr);
      }
      /* Numerical flux of By at cell corner along Y to calculate Ez */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	int s1l,s1r;
	int sl,sr;
	int sm2,sm,sp,sp2;
	double *ul,*ur;
	ul=(double*)malloc(sizeof(double)*2*ny);
	ur=(double*)malloc(sizeof(double)*2*ny);
	for (i=3;i<nx-2;i++){
	  /* Interpolation */
#pragma simd
	  for (j=2;j<ny-2;j++){
	    s1l=j+1;
	    s1r=j;
	    ss=nx*(ny*k+j)+i;
	    sl=nx*(ny*k+(j+1))+i;
	    sr=nx*(ny*k+j)+i;
	    sm2=nx*(ny*k+(j-2))+i;
	    sm=nx*(ny*k+(j-1))+i;
	    sp=nx*(ny*k+(j+1))+i;
	    sp2=nx*(ny*k+(j+2))+i;
	    double val[]={fx[8*sm2+5],fx[8*sm+5],fx[8*ss+5],fx[8*sp+5],fx[8*sp2+5]};
	    double vac[]={fc[2*sm2+0],fc[2*sm+0],fc[2*ss+0],fc[2*sp+0],fc[2*sp2+0]};

	    double vl[2],vr[2];
	    lfun_lr(val,&vl[0],&vr[0]); /* Upwind part */
	    lfun_lr(vac,&vl[1],&vr[1]); /* Central part */
	    ul[2*s1l+0]=vl[0];
	    ur[2*s1r+0]=vr[0];
	    ul[2*s1l+1]=vl[1];
	    ur[2*s1r+1]=vr[1];
	  }
	  /* CT method */
#pragma simd
	  for (j=3;j<ny-2;j++){
	    s1=j;
	    ss=nx*(ny*k+j)+i;
	    ez[ss]+=-0.5*((ul[2*s1+0]+ur[2*s1+0])+(1.0-ct[nxyz*0+ss])*(ul[2*s1+1]+ur[2*s1+1]));
	  }
	}
	free(ul);
	free(ur);
      }
      /* Numerical flux of Bz at cell corner along Z to calculate Ey */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
	int s1l,s1r;
	int sl,sr;
	int sm2,sm,sp,sp2;
	double *ul,*ur;
	ul=(double*)malloc(sizeof(double)*2*nz);
	ur=(double*)malloc(sizeof(double)*2*nz);
	for (i=3;i<nx-2;i++){
	  /* Interpolation */
#pragma simd
	  for (k=2;k<nz-2;k++){
	    s1l=k+1;
	    s1r=k;
	    ss=nx*(ny*k+j)+i;
	    sl=nx*(ny*(k+1)+j)+i;
	    sr=nx*(ny*k+j)+i;
	    sm2=nx*(ny*(k-2)+j)+i;
	    sm=nx*(ny*(k-1)+j)+i;
	    sp=nx*(ny*(k+1)+j)+i;
	    sp2=nx*(ny*(k+2)+j)+i;
	    double val[]={fx[8*sm2+6],fx[8*sm+6],fx[8*ss+6],fx[8*sp+6],fx[8*sp2+6]};
	    double vac[]={fc[2*sm2+1],fc[2*sm+1],fc[2*ss+1],fc[2*sp+1],fc[2*sp2+1]};

	    double vl[2],vr[2];
	    lfun_lr(val,&vl[0],&vr[0]); /* Upwind part */
	    lfun_lr(vac,&vl[1],&vr[1]); /* Central part */
	    ul[2*s1l+0]=vl[0];
	    ur[2*s1r+0]=vr[0];
	    ul[2*s1l+1]=vl[1];
	    ur[2*s1r+1]=vr[1];
	  }
	  /* CT method */
#pragma simd
	  for (k=3;k<nz-2;k++){
	    s1=k;
	    ss=nx*(ny*k+j)+i;
	    ey[ss]+=+0.5*((ul[2*s1+0]+ur[2*s1+0])+ct[nxyz*2+ss]*(ul[2*s1+1]+ur[2*s1+1]));
	  }
	}
	free(ul);
	free(ur);
      }

      /* Calculate numerical flux at cell face along Y */	  
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	int s1l,s1r;
	int sl,sr;
	int sm2,sm,sp,sp2;
	double *ul,*ur,*ql,*qr;
	ul=(double*)malloc(sizeof(double)*8*ny);
	ur=(double*)malloc(sizeof(double)*8*ny);
	ql=(double*)malloc(sizeof(double)*2*ny);
	qr=(double*)malloc(sizeof(double)*2*ny);
	for (i=0;i<nx;i++){
	  /* Interpolation */
	  for (j=2;j<ny-2;j++){
	    s1l=j+1;
	    s1r=j;
	    ss=nx*(ny*k+j)+i;
	    sl=nx*(ny*k+(j+1))+i;
	    sr=nx*(ny*k+j)+i;
	    sm2=nx*(ny*k+(j-2))+i;
	    sm=nx*(ny*k+(j-1))+i;
	    sp=nx*(ny*k+(j+1))+i;
	    sp2=nx*(ny*k+(j+2))+i;
	    /* Primitive variables in the stencil */
	    double ros[]={ro[sm2],ro[sm],ro[ss],ro[sp],ro[sp2]};
	    double vxs[]={vx[sm2],vx[sm],vx[ss],vx[sp],vx[sp2]};
	    double vys[]={vy[sm2],vy[sm],vy[ss],vy[sp],vy[sp2]};
	    double vzs[]={vz[sm2],vz[sm],vz[ss],vz[sp],vz[sp2]};
	    double bxs[]={cx[sm2],cx[sm],cx[ss],cx[sp],cx[sp2]};
	    double bys[]={cy[sm2],cy[sm],cy[ss],cy[sp],cy[sp2]};
	    double bzs[]={cz[sm2],cz[sm],cz[ss],cz[sp],cz[sp2]};
	    double prs[]={pr[sm2],pr[sm],pr[ss],pr[sp],pr[sp2]};

	    double vl[7],vr[7];

	    if (max(m2[sm],max(m2[ss],m2[sp])) >= 1.0){
	      mhd_c_reconst(ros,vys,vzs,vxs,bzs,bxs,prs,bys[ns/2],gamma,ns,vl,vr,func_lr);
	    } else{
	      mhd_a_reconst(ros,vys,vzs,vxs,bzs,bxs,prs,bys[ns/2],gamma,ns,vl,vr,func_lr);
	    }

	    /* Left-face @ j+1/2 */
	    ul[8*s1l+0]=vl[0];	/* ro */
	    ul[8*s1l+1]=vl[1];	/* vy */
	    ul[8*s1l+2]=vl[2];	/* vz */
	    ul[8*s1l+3]=vl[3];	/* vx */
	    ul[8*s1l+4]=by[sl];	/* by */
	    ul[8*s1l+5]=vl[4];	/* bz */
	    ul[8*s1l+6]=vl[5];	/* bx */
	    ul[8*s1l+7]=vl[6];	/* pr */
	    /* Right-face @ j-1/2 */
	    ur[8*s1r+0]=vr[0];	/* ro */
	    ur[8*s1r+1]=vr[1];	/* vy */
	    ur[8*s1r+2]=vr[2];	/* vz */
	    ur[8*s1r+3]=vr[3];	/* vx */
	    ur[8*s1r+4]=by[sr];	/* by */
	    ur[8*s1r+5]=vr[4];	/* bz */
	    ur[8*s1r+6]=vr[5];	/* bx */
	    ur[8*s1r+7]=vr[6];	/* pr */

	    /* Linear interpolation of numerical flux of Bz */
	    double vaz[]={cz[sm2]*vy[sm2],cz[sm]*vy[sm],cz[ss]*vy[ss],cz[sp]*vy[sp],cz[sp2]*vy[sp2]};
	    lfun_lr(vaz,&vl[0],&vr[0]);
	    lfun_lr(vzs,&vl[1],&vr[1]);
	    ql[2*s1l+0]=vl[0]-by[sl]*vl[1];
	    qr[2*s1r+0]=vr[0]-by[sr]*vr[1];

	    /* Linear interpolation of numerical flux of Bx */
	    double vax[]={cx[sm2]*vy[sm2],cx[sm]*vy[sm],cx[ss]*vy[ss],cx[sp]*vy[sp],cx[sp2]*vy[sp2]};
	    lfun_lr(vax,&vl[0],&vr[0]);
	    lfun_lr(vxs,&vl[1],&vr[1]);
	    ql[2*s1l+1]=vl[0]-by[sl]*vl[1];
	    qr[2*s1r+1]=vr[0]-by[sr]*vr[1];
	  }
	  /* Riemann solver */
#pragma simd
	  for (j=3;j<ny-2;j++){
	    double bn,flux[8];
	    s1=j;
	    ss=nx*(ny*k+j)+i;
	    bn=0.5*(ul[8*s1+4]+ur[8*s1+4]);

	    double dvtm=min(min(dvz[nx*(ny*k+(j-1))+i],dvz[nx*(ny*k+j)+i]),
			    min(dvx[nx*(ny*k+(j-1))+i],dvx[nx*(ny*k+j)+i]));
	    double dvsd[2]={(vy[nx*(ny*k+j)+i]-vy[nx*(ny*k+(j-1))+i]),dvtm};
	    func_flux(ul[8*s1+0],ul[8*s1+1],ul[8*s1+2],ul[8*s1+3],ul[8*s1+5],ul[8*s1+6],ul[8*s1+7],
		      ur[8*s1+0],ur[8*s1+1],ur[8*s1+2],ur[8*s1+3],ur[8*s1+5],ur[8*s1+6],ur[8*s1+7],
		      bn,gamma,dvsd,
		      &flux[0],&flux[1],&flux[2],&flux[3],&flux[5],&flux[6],&flux[7]);

	    fy[8*ss+0]=flux[0];	/* ro */
	    fy[8*ss+2]=flux[1];	/* my */
	    fy[8*ss+3]=flux[2];	/* mz */
	    fy[8*ss+1]=flux[3];	/* mx */
	    fy[8*ss+5]=0;	/* by */
	    fy[8*ss+6]=flux[5];	/* bz */
	    fy[8*ss+4]=flux[6];	/* bx */
	    fy[8*ss+7]=flux[7];	/* en */

	    /* Divide central and upwind parts in numerical flux of Bz */
	    fc[2*ss+0]=0.5*(ql[2*s1+0]+qr[2*s1+0]); /* Central part */
	    fy[8*ss+6]-=fc[2*ss+0];			     /* Upwind part */

	    /* Divide central and upwind parts in numerical flux of Bx */
	    fc[2*ss+1]=0.5*(ql[2*s1+1]+qr[2*s1+1]); /* Central part */
	    fy[8*ss+4]-=fc[2*ss+1];			     /* Upwind part */
	  }
	}
	free(ul);
	free(ur);
	free(ql);
	free(qr);
      }
      /* Numerical flux of Bz at cell corner along Z to calculate Ex */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
	int s1l,s1r;
	int sl,sr;
	int sm2,sm,sp,sp2;
	double *ul,*ur;
	ul=(double*)malloc(sizeof(double)*2*nz);
	ur=(double*)malloc(sizeof(double)*2*nz);
	for (i=0;i<nx;i++){
	  /* Interpolation */
#pragma simd
	  for (k=2;k<nz-2;k++){
	    s1l=k+1;
	    s1r=k;
	    ss=nx*(ny*k+j)+i;
	    sl=nx*(ny*(k+1)+j)+i;
	    sr=nx*(ny*k+j)+i;
	    sm2=nx*(ny*(k-2)+j)+i;
	    sm=nx*(ny*(k-1)+j)+i;
	    sp=nx*(ny*(k+1)+j)+i;
	    sp2=nx*(ny*(k+2)+j)+i;
	    double val[]={fy[8*sm2+6],fy[8*sm+6],fy[8*ss+6],fy[8*sp+6],fy[8*sp2+6]};
	    double vac[]={fc[2*sm2+0],fc[2*sm+0],fc[2*ss+0],fc[2*sp+0],fc[2*sp2+0]};

	    double vl[2],vr[2];
	    lfun_lr(val,&vl[0],&vr[0]); /* Upwind part */
	    lfun_lr(vac,&vl[1],&vr[1]); /* Central part */
	    ul[2*s1l+0]=vl[0];
	    ur[2*s1r+0]=vr[0];
	    ul[2*s1l+1]=vl[1];
	    ur[2*s1r+1]=vr[1];
	  }
	  /* CT method */
#pragma simd
	  for (k=3;k<nz-2;k++){
	    s1=k;
	    ss=nx*(ny*k+j)+i;
	    ex[ss]+=-0.5*((ul[2*s1+0]+ur[2*s1+0])+(1.0-ct[nxyz*1+ss])*(ul[2*s1+1]+ur[2*s1+1]));
	  }
	}
	free(ul);
	free(ur);
      }
      /* Numerical flux of Bx at cell corner along X to calculate Ez */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	int s1l,s1r;
	int sl,sr;
	int sm2,sm,sp,sp2;
	double *ul,*ur;
	ul=(double*)malloc(sizeof(double)*2*nx);
	ur=(double*)malloc(sizeof(double)*2*nx);
	for (j=3;j<ny-2;j++){
	  /* Interpolation */
#pragma simd
	  for (i=2;i<nx-2;i++){
	    s1l=i+1;
	    s1r=i;
	    ss=nx*(ny*k+j)+i;
	    sl=nx*(ny*k+j)+(i+1);
	    sr=nx*(ny*k+j)+i;
	    sm2=nx*(ny*k+j)+(i-2);
	    sm=nx*(ny*k+j)+(i-1);
	    sp=nx*(ny*k+j)+(i+1);
	    sp2=nx*(ny*k+j)+(i+2);
	    double val[]={fy[8*sm2+4],fy[8*sm+4],fy[8*ss+4],fy[8*sp+4],fy[8*sp2+4]};
	    double vac[]={fc[2*sm2+1],fc[2*sm+1],fc[2*ss+1],fc[2*sp+1],fc[2*sp2+1]};

	    double vl[2],vr[2];
	    lfun_lr(val,&vl[0],&vr[0]); /* Upwind part */
	    lfun_lr(vac,&vl[1],&vr[1]); /* Central part */
	    ul[2*s1l+0]=vl[0];
	    ur[2*s1r+0]=vr[0];
	    ul[2*s1l+1]=vl[1];
	    ur[2*s1r+1]=vr[1];
	  }
	  /* CT method */
#pragma simd
	  for (i=3;i<nx-2;i++){
	    s1=i;
	    ss=nx*(ny*k+j)+i;
	    ez[ss]+=+0.5*((ul[2*s1+0]+ur[2*s1+0])+ct[nxyz*0+ss]*(ul[2*s1+1]+ur[2*s1+1]));
	  }
	}
	free(ul);
	free(ur);
      }

      /* Calculate numerical flux at cell face along Z */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
	int s1l,s1r;
	int sl,sr;
	int sm2,sm,sp,sp2;
	double *ul,*ur,*ql,*qr;
	ul=(double*)malloc(sizeof(double)*8*nz);
	ur=(double*)malloc(sizeof(double)*8*nz);
	ql=(double*)malloc(sizeof(double)*2*nz);
	qr=(double*)malloc(sizeof(double)*2*nz);
	for (i=0;i<nx;i++){
	  /* Interpolation */
	  for (k=2;k<nz-2;k++){
	    s1l=k+1;
	    s1r=k;
	    ss=nx*(ny*k+j)+i;
	    sl=nx*(ny*(k+1)+j)+i;
	    sr=nx*(ny*k+j)+i;
	    sm2=nx*(ny*(k-2)+j)+i;
	    sm=nx*(ny*(k-1)+j)+i;
	    sp=nx*(ny*(k+1)+j)+i;
	    sp2=nx*(ny*(k+2)+j)+i;
	    /* Primitive variables in the stencil */
	    double ros[]={ro[sm2],ro[sm],ro[ss],ro[sp],ro[sp2]};
	    double vxs[]={vx[sm2],vx[sm],vx[ss],vx[sp],vx[sp2]};
	    double vys[]={vy[sm2],vy[sm],vy[ss],vy[sp],vy[sp2]};
	    double vzs[]={vz[sm2],vz[sm],vz[ss],vz[sp],vz[sp2]};
	    double bxs[]={cx[sm2],cx[sm],cx[ss],cx[sp],cx[sp2]};
	    double bys[]={cy[sm2],cy[sm],cy[ss],cy[sp],cy[sp2]};
	    double bzs[]={cz[sm2],cz[sm],cz[ss],cz[sp],cz[sp2]};
	    double prs[]={pr[sm2],pr[sm],pr[ss],pr[sp],pr[sp2]};

	    double vl[7],vr[7];

	    if (max(m2[sm],max(m2[ss],m2[sp])) >= 1.0){
	      mhd_c_reconst(ros,vzs,vxs,vys,bxs,bys,prs,bzs[ns/2],gamma,ns,vl,vr,func_lr);
	    } else{
	      mhd_a_reconst(ros,vzs,vxs,vys,bxs,bys,prs,bzs[ns/2],gamma,ns,vl,vr,func_lr);
	    }

	    /* Left-face @ k+1/2 */
	    ul[8*s1l+0]=vl[0];	/* ro */
	    ul[8*s1l+1]=vl[1];	/* vz */
	    ul[8*s1l+2]=vl[2];	/* vx */
	    ul[8*s1l+3]=vl[3];	/* vy */
	    ul[8*s1l+4]=bz[sl];	/* bz */
	    ul[8*s1l+5]=vl[4];	/* bx */
	    ul[8*s1l+6]=vl[5];	/* by */
	    ul[8*s1l+7]=vl[6];	/* pr */
	    /* Right-face @ k-1/2 */
	    ur[8*s1r+0]=vr[0];	/* ro */
	    ur[8*s1r+1]=vr[1];	/* vz */
	    ur[8*s1r+2]=vr[2];	/* vx */
	    ur[8*s1r+3]=vr[3];	/* vy */
	    ur[8*s1r+4]=bz[sr];	/* bz */
	    ur[8*s1r+5]=vr[4];	/* bx */
	    ur[8*s1r+6]=vr[5];	/* by */
	    ur[8*s1r+7]=vr[6];	/* pr */

	    /* Linear interpolation of numerical flux of Bx */
	    double vax[]={cx[sm2]*vz[sm2],cx[sm]*vz[sm],cx[ss]*vz[ss],cx[sp]*vz[sp],cx[sp2]*vz[sp2]};
	    lfun_lr(vax,&vl[0],&vr[0]);
	    lfun_lr(vxs,&vl[1],&vr[1]);
	    ql[2*s1l+0]=vl[0]-bz[sl]*vl[1];
	    qr[2*s1r+0]=vr[0]-bz[sr]*vr[1];

	    /* Linear interpolation of numerical flux of By */
	    double vay[]={cy[sm2]*vz[sm2],cy[sm]*vz[sm],cy[ss]*vz[ss],cy[sp]*vz[sp],cy[sp2]*vz[sp2]};
	    lfun_lr(vay,&vl[0],&vr[0]);
	    lfun_lr(vys,&vl[1],&vr[1]);
	    ql[2*s1l+1]=vl[0]-bz[sl]*vl[1];
	    qr[2*s1r+1]=vr[0]-bz[sr]*vr[1];
	  }
	  /* Riemann solver */
#pragma simd
	  for (k=3;k<nz-2;k++){
	    double bn,flux[8];
	    s1=k;
	    ss=nx*(ny*k+j)+i;
	    bn=0.5*(ul[8*s1+4]+ur[8*s1+4]);

	    double dvtm=min(min(dvx[nx*(ny*(k-1)+j)+i],dvx[nx*(ny*k+j)+i]),
			    min(dvy[nx*(ny*(k-1)+j)+i],dvy[nx*(ny*k+j)+i]));
	    double dvsd[2]={(vz[nx*(ny*k+j)+i]-vz[nx*(ny*(k-1)+j)+i]),dvtm};
	    func_flux(ul[8*s1+0],ul[8*s1+1],ul[8*s1+2],ul[8*s1+3],ul[8*s1+5],ul[8*s1+6],ul[8*s1+7],
		      ur[8*s1+0],ur[8*s1+1],ur[8*s1+2],ur[8*s1+3],ur[8*s1+5],ur[8*s1+6],ur[8*s1+7],
		      bn,gamma,dvsd,
		      &flux[0],&flux[1],&flux[2],&flux[3],&flux[5],&flux[6],&flux[7]);

	    fz[8*ss+0]=flux[0];	/* ro */
	    fz[8*ss+3]=flux[1];	/* mz */
	    fz[8*ss+1]=flux[2];	/* mx */
	    fz[8*ss+2]=flux[3];	/* my */
	    fz[8*ss+6]=0;	/* bz */
	    fz[8*ss+4]=flux[5];	/* bx */
	    fz[8*ss+5]=flux[6];	/* by */
	    fz[8*ss+7]=flux[7];	/* en */

	    /* Divide central and upwind parts in numerical flux of Bx */
	    fc[2*ss+0]=0.5*(ql[2*s1+0]+qr[2*s1+0]); /* Central part */
	    fz[8*ss+4]-=fc[2*ss+0];			     /* Upwind part */

	    /* Divide central and upwind parts in numerical flux of By */
	    fc[2*ss+1]=0.5*(ql[2*s1+1]+qr[2*s1+1]); /* Central part */
	    fz[8*ss+5]-=fc[2*ss+1];			     /* Upwind part */
	  }
	}
	free(ul);
	free(ur);
	free(ql);
	free(qr);
      }
      /* Numerical flux of Bx at cell corner along X to calculate Ey */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=3;k<nz-2;k++){
	int s1l,s1r;
	int sl,sr;
	int sm2,sm,sp,sp2;
	double *ul,*ur;
	ul=(double*)malloc(sizeof(double)*2*nx);
	ur=(double*)malloc(sizeof(double)*2*nx);
	for (j=0;j<ny;j++){
	  /* Interpolation */
#pragma simd
	  for (i=2;i<nx-2;i++){
	    s1l=i+1;
	    s1r=i;
	    ss=nx*(ny*k+j)+i;
	    sl=nx*(ny*k+j)+(i+1);
	    sr=nx*(ny*k+j)+i;
	    sm2=nx*(ny*k+j)+(i-2);
	    sm=nx*(ny*k+j)+(i-1);
	    sp=nx*(ny*k+j)+(i+1);
	    sp2=nx*(ny*k+j)+(i+2);
	    double val[]={fz[8*sm2+4],fz[8*sm+4],fz[8*ss+4],fz[8*sp+4],fz[8*sp2+4]};
	    double vac[]={fc[2*sm2+0],fc[2*sm+0],fc[2*ss+0],fc[2*sp+0],fc[2*sp2+0]};

	    double vl[2],vr[2];
	    lfun_lr(val,&vl[0],&vr[0]); /* Upwind part */
	    lfun_lr(vac,&vl[1],&vr[1]); /* Central part */
	    ul[2*s1l+0]=vl[0];
	    ur[2*s1r+0]=vr[0];
	    ul[2*s1l+1]=vl[1];
	    ur[2*s1r+1]=vr[1];
	  }
	  /* CT method */
#pragma simd
	  for (i=3;i<nx-2;i++){
	    s1=i;
	    ss=nx*(ny*k+j)+i;
	    ey[ss]+=-0.5*((ul[2*s1+0]+ur[2*s1+0])+(1.0-ct[nxyz*2+ss])*(ul[2*s1+1]+ur[2*s1+1]));
	  }
	}
	free(ul);
	free(ur);
      }
      /* Numerical flux of By at cell corner along Y to calculate Ex */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=3;k<nz-2;k++){
	int s1l,s1r;
	int sl,sr;
	int sm2,sm,sp,sp2;
	double *ul,*ur;
	ul=(double*)malloc(sizeof(double)*2*ny);
	ur=(double*)malloc(sizeof(double)*2*ny);
	for (i=0;i<nx;i++){
	  /* Interpolation */
#pragma simd
	  for (j=2;j<ny-2;j++){
	    s1l=j+1;
	    s1r=j;
	    ss=nx*(ny*k+j)+i;
	    sl=nx*(ny*k+(j+1))+i;
	    sr=nx*(ny*k+j)+i;
	    sm2=nx*(ny*k+(j-2))+i;
	    sm=nx*(ny*k+(j-1))+i;
	    sp=nx*(ny*k+(j+1))+i;
	    sp2=nx*(ny*k+(j+2))+i;
	    double val[]={fz[8*sm2+5],fz[8*sm+5],fz[8*ss+5],fz[8*sp+5],fz[8*sp2+5]};
	    double vac[]={fc[2*sm2+1],fc[2*sm+1],fc[2*ss+1],fc[2*sp+1],fc[2*sp2+1]};

	    double vl[2],vr[2];
	    lfun_lr(val,&vl[0],&vr[0]); /* Upwind part */
	    lfun_lr(vac,&vl[1],&vr[1]); /* Central part */
	    ul[2*s1l+0]=vl[0];
	    ur[2*s1r+0]=vr[0];
	    ul[2*s1l+1]=vl[1];
	    ur[2*s1r+1]=vr[1];
	  }
	  /* CT method */
#pragma simd
	  for (j=3;j<ny-2;j++){
	    s1=j;
	    ss=nx*(ny*k+j)+i;
	    ex[ss]+=+0.5*((ul[2*s1+0]+ur[2*s1+0])+ct[nxyz*1+ss]*(ul[2*s1+1]+ur[2*s1+1]));
	  }
	}
	free(ul);
	free(ur);
      }

      /* Update variable at cell center */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (k=zoff;k<nz-zoff;k++){
	int sip,sjp,skp,si0,sj0,sk0;
	int sim,sjm,skm,sip2,sjp2,skp2;
	for (j=yoff;j<ny-yoff;j++){
#pragma simd
	  for (i=xoff;i<nx-xoff;i++){
	    ss=nx*(ny*k+j)+i;
	    sip=nx*(ny*k+j)+(i+1);
	    sjp=nx*(ny*k+(j+1))+i;
	    skp=nx*(ny*(k+1)+j)+i;
	    si0=nx*(ny*k+j)+i;
	    sj0=nx*(ny*k+j)+i;
	    sk0=nx*(ny*k+j)+i;
	    sim=nx*(ny*k+j)+(i-1);
	    sjm=nx*(ny*k+(j-1))+i;
	    skm=nx*(ny*(k-1)+j)+i;
	    sip2=nx*(ny*k+j)+(i+2);
	    sjp2=nx*(ny*k+(j+2))+i;
	    skp2=nx*(ny*(k+2)+j)+i;
	    double du[3][8];

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
	    du[1][7]=func_df(&valy[7][0]);
	    /* dH/dz */
	    double valz[8][4]={{fz[8*skm+0],fz[8*sk0+0],fz[8*skp+0],fz[8*skp2+0]},
			       {fz[8*skm+1],fz[8*sk0+1],fz[8*skp+1],fz[8*skp2+1]},
			       {fz[8*skm+2],fz[8*sk0+2],fz[8*skp+2],fz[8*skp2+2]},
			       {fz[8*skm+3],fz[8*sk0+3],fz[8*skp+3],fz[8*skp2+3]},
			       {fz[8*skm+4],fz[8*sk0+4],fz[8*skp+4],fz[8*skp2+4]},
			       {fz[8*skm+5],fz[8*sk0+5],fz[8*skp+5],fz[8*skp2+5]},
			       {fz[8*skm+6],fz[8*sk0+6],fz[8*skp+6],fz[8*skp2+6]},
			       {fz[8*skm+7],fz[8*sk0+7],fz[8*skp+7],fz[8*skp2+7]}};
	    du[2][0]=func_df(&valz[0][0]);
	    du[2][1]=func_df(&valz[1][0]);
	    du[2][2]=func_df(&valz[2][0]);
	    du[2][3]=func_df(&valz[3][0]);
	    du[2][7]=func_df(&valz[7][0]);
	  
	    ro[ss]=rk_fac[rk][0]*ut[8*ss+0]+rk_fac[rk][1]*(ro[ss]-dtdx*du[0][0]-dtdy*du[1][0]-dtdz*du[2][0]);
	    mx[ss]=rk_fac[rk][0]*ut[8*ss+1]+rk_fac[rk][1]*(mx[ss]-dtdx*du[0][1]-dtdy*du[1][1]-dtdz*du[2][1]);
	    my[ss]=rk_fac[rk][0]*ut[8*ss+2]+rk_fac[rk][1]*(my[ss]-dtdx*du[0][2]-dtdy*du[1][2]-dtdz*du[2][2]);
	    mz[ss]=rk_fac[rk][0]*ut[8*ss+3]+rk_fac[rk][1]*(mz[ss]-dtdx*du[0][3]-dtdy*du[1][3]-dtdz*du[2][3]);
	    en[ss]=rk_fac[rk][0]*ut[8*ss+7]+rk_fac[rk][1]*(en[ss]-dtdx*du[0][7]-dtdy*du[1][7]-dtdz*du[2][7]);
	  }
	}
      }
      /* Update CT Bx and By */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (k=zoff;k<nz-zoff;k++){
	int sip,sjp,skp,si0,sj0,sk0;
	int sim,sjm,skm,sip2,sjp2,skp2;
	for (j=yoff;j<ny-yoff;j++){
#pragma simd
	  for (i=xoff;i<nx-xoff+1;i++){
	    ss=nx*(ny*k+j)+i;
	    sjp=nx*(ny*k+(j+1))+i;
	    skp=nx*(ny*(k+1)+j)+i;
	    sj0=nx*(ny*k+j)+i;
	    sk0=nx*(ny*k+j)+i;
	    sjm=nx*(ny*k+(j-1))+i;
	    skm=nx*(ny*(k-1)+j)+i;
	    sjp2=nx*(ny*k+(j+2))+i;
	    skp2=nx*(ny*(k+2)+j)+i;
	    double du[2];
	    
	    double val[2][4]={{ez[sjm],ez[sj0],ez[sjp],ez[sjp2]},
			      {ey[skm],ey[sk0],ey[skp],ey[skp2]}};
	    du[0]=func_df(&val[0][0]);
	    du[1]=func_df(&val[1][0]);
	    bx[ss]=rk_fac[rk][0]*ut[8*ss+4]+rk_fac[rk][1]*(bx[ss]-dtdy*du[0]+dtdz*du[1]);
	  }
	}
	for (j=yoff;j<ny-yoff+1;j++){
#pragma simd
	  for (i=xoff;i<nx-xoff;i++){
	    ss=nx*(ny*k+j)+i;
	    sip=nx*(ny*k+j)+(i+1);
	    skp=nx*(ny*(k+1)+j)+i;
	    si0=nx*(ny*k+j)+i;
	    sk0=nx*(ny*k+j)+i;
	    sim=nx*(ny*k+j)+(i-1);
	    skm=nx*(ny*(k-1)+j)+i;
	    sip2=nx*(ny*k+j)+(i+2);
	    skp2=nx*(ny*(k+2)+j)+i;
	    double du[2];
	    
	    double val[2][4]={{ex[skm],ex[sk0],ex[skp],ex[skp2]},
			      {ez[sim],ez[si0],ez[sip],ez[sip2]}};
	    du[0]=func_df(&val[0][0]);
	    du[1]=func_df(&val[1][0]);
	    by[ss]=rk_fac[rk][0]*ut[8*ss+5]+rk_fac[rk][1]*(by[ss]-dtdz*du[0]+dtdx*du[1]);
	  }
	}
      }
      /* Update CT Bz */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=zoff;k<nz-zoff+1;k++){
	int sip,sjp,skp,si0,sj0,sk0;
	int sim,sjm,skm,sip2,sjp2,skp2;
	for (j=yoff;j<ny-yoff;j++){
#pragma simd
	  for (i=xoff;i<nx-xoff;i++){
	    ss=nx*(ny*k+j)+i;
	    sip=nx*(ny*k+j)+(i+1);
	    sjp=nx*(ny*k+(j+1))+i;
	    si0=nx*(ny*k+j)+i;
	    sj0=nx*(ny*k+j)+i;
	    sim=nx*(ny*k+j)+(i-1);
	    sjm=nx*(ny*k+(j-1))+i;
	    sip2=nx*(ny*k+j)+(i+2);
	    sjp2=nx*(ny*k+(j+2))+i;
	    double du[2];

	    double val[2][4]={{ey[sim],ey[si0],ey[sip],ey[sip2]},
			      {ex[sjm],ex[sj0],ex[sjp],ex[sjp2]}};
	    du[0]=func_df(&val[0][0]);
	    du[1]=func_df(&val[1][0]);
	    bz[ss]=rk_fac[rk][0]*ut[8*ss+6]+rk_fac[rk][1]*(bz[ss]-dtdx*du[0]+dtdy*du[1]);
	  }
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
    mpi_sdrv3d(p,8,nx,ny,nz,xoff,yoff,zoff,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_xbc3d(&ro[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_xbc3d(&mx[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_xbc3d(&my[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_xbc3d(&mz[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_xbc3d(&en[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_xbc3d(&bx[0],nx,ny,nz,xoff,yoff,zoff,1,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_xbc3d(&by[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_xbc3d(&bz[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_ybc3d(&ro[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_ybc3d(&mx[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_ybc3d(&my[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_ybc3d(&mz[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_ybc3d(&en[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_ybc3d(&bx[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_ybc3d(&by[0],nx,ny,nz,xoff,yoff,zoff,1,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_ybc3d(&bz[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_zbc3d(&ro[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_zbc3d(&mx[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_zbc3d(&my[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_zbc3d(&mz[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_zbc3d(&en[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_zbc3d(&bx[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_zbc3d(&by[0],nx,ny,nz,xoff,yoff,zoff,0,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_zbc3d(&bz[0],nx,ny,nz,xoff,yoff,zoff,1,+0,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
  }
  
  free(ut);
  free(fx);
  free(fy);
  free(fz);
  free(vx);
  free(vy);
  free(vz);
  free(pr);
  free(cx);
  free(cy);
  free(cz);
  free(ex);
  free(ey);
  free(ez);
  free(ct);
  free(fc);
  free(m2);
  free(dvx);
  free(dvy);
  free(dvz);
}
