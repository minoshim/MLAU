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

#ifndef _MHD_EIGEN_H_
#define _MHD_EIGEN_H_

#include <math.h>
#include "funcs.h"

inline void mhd_eigen_st(double ro, double vx, double vy, double vz, double bx, double by, double bz, double pr, double gamma,
			 double *al1, double *al2, double *al3, double *al4, double *al5, double *al6, double *al7,
			 double *l1, double *l2, double *l3, double *l4, double *l5, double *l6, double *l7, 
			 double *r1, double *r2, double *r3, double *r4, double *r5, double *r6, double *r7)
/* Calculate eigenvalues and eigenvectors for non-conservative form, based on ATHENA */
/* Stone+ 08 ApJS, 178,137 */
/* Outputs: */
/* al?: eigenvalues (NOT absolute) */
/* l?: left eigenvectors (row vectors) */
/* r?: right eigenvectors (column vectors) */
{
  double eps=1e-12;
  double sqrt2=sqrt(2.0);
  double isqrt2=1.0/sqrt2;

  /* Wave speeds */
  double iro=1.0/ro;
  double sro=sqrt(ro);
  double isro=1.0/sro;
  double bxsr=bx*isro;
  double bysr=by*isro;
  double bzsr=bz*isro;
  double bx2;
  double bt2;
  double b2;
  bx2=bxsr*bxsr;
  bt2=bysr*bysr+bzsr*bzsr;
  b2=bx2+bt2;
  double a2=fabs(gamma*pr/ro);
  double aa=sqrt(a2);
  double ca2=bx2;
  double ca=sqrt(ca2);
  double junk;
  junk=sqrt(fabs((a2+b2)*(a2+b2)-4.0*a2*bx2));
  double cf2=fabs(0.5*(a2+b2+junk));
  double cf=sqrt(cf2);
  double cs2=fabs(0.5*(a2+b2-junk));
  double cs=sqrt(cs2);

  /* Quantities for engenvector calculation */
  bx2=bx*bx;
  bt2=by*by+bz*bz;
  b2=bx2+bt2;
  double sbt2=sqrt(bt2);
  double isbt2=1.0/sbt2;
  int flag;
  flag=(bt2 > eps*b2);
  double betay=(flag)?(by*isbt2):(0.5*sqrt(2.0));
  double betaz=(flag)?(bz*isbt2):(0.5*sqrt(2.0));
  double scfcs=sqrt(fabs(cf2-cs2));
  double iscfcs=1.0/scfcs;
  flag=((bt2 > eps*b2) || (fabs(gamma*pr-bx2) > eps*gamma*pr));
  double alphaf=(flag)?(sqrt(fabs(a2-cs2))*iscfcs):1.0;
  double alphas=(flag)?(sqrt(fabs(cf2-a2))*iscfcs):0.0;
  double sgnbx=(bx > 0)?(1.0):(-1.0);
  double cfalf=cf*alphaf;
  double csals=cs*alphas;
  double aaalf=aa*alphaf;
  double aaals=aa*alphas;
  double denom=0.5/a2;

  /* Eigenvalues */
  *al1=(vx-cf);
  *al2=(vx-ca);
  *al3=(vx-cs);
  *al4=(vx);
  *al5=(vx+cs);
  *al6=(vx+ca);
  *al7=(vx+cf);

  /* Left eigenvectors */
  l1[0]=0;
  l1[1]=-cfalf*denom;
  l1[2]=+csals*denom*sgnbx*betay;
  l1[3]=+csals*denom*sgnbx*betaz;
  l1[4]=+aaals*denom*isro*betay;
  l1[5]=+aaals*denom*isro*betaz;
  l1[6]=alphaf*denom*iro;
  l7[0]=+l1[0];
  l7[1]=-l1[1];
  l7[2]=-l1[2];
  l7[3]=-l1[3];
  l7[4]=+l1[4];
  l7[5]=+l1[5];
  l7[6]=+l1[6];
  l2[0]=0;
  l2[1]=0;
  l2[2]=-0.5*betaz;
  l2[3]=+0.5*betay;
  l2[4]=-0.5*isro*betaz*sgnbx;
  l2[5]=+0.5*isro*betay*sgnbx;
  l2[6]=0;
  l6[0]=+l2[0];
  l6[1]=+l2[1];
  l6[2]=-l2[2];
  l6[3]=-l2[3];
  l6[4]=+l2[4];
  l6[5]=+l2[5];
  l6[6]=+l2[6];
  l3[0]=0;
  l3[1]=-csals*denom;
  l3[2]=-cfalf*denom*sgnbx*betay;
  l3[3]=-cfalf*denom*sgnbx*betaz;
  l3[4]=-aaalf*denom*isro*betay;
  l3[5]=-aaalf*denom*isro*betaz;
  l3[6]=alphas*denom*iro;
  l5[0]=+l3[0];
  l5[1]=-l3[1];
  l5[2]=-l3[2];
  l5[3]=-l3[3];
  l5[4]=+l3[4];
  l5[5]=+l3[5];
  l5[6]=+l3[6];
  l4[0]=1.0;
  l4[1]=0;
  l4[2]=0;
  l4[3]=0;
  l4[4]=0;
  l4[5]=0;
  l4[6]=-1.0/a2;

  /* Right eigenvectors */
  r1[0]=ro*alphaf;
  r1[1]=-cfalf;
  r1[2]=+csals*sgnbx*betay;
  r1[3]=+csals*sgnbx*betaz;
  r1[4]=+aaals*sro*betay;
  r1[5]=+aaals*sro*betaz;
  r1[6]=alphaf*gamma*pr;
  r7[0]=+r1[0];
  r7[1]=-r1[1];
  r7[2]=-r1[2];
  r7[3]=-r1[3];
  r7[4]=+r1[4];
  r7[5]=+r1[5];
  r7[6]=+r1[6];
  r2[0]=0;
  r2[1]=0;
  r2[2]=-betaz;
  r2[3]=+betay;
  r2[4]=-sro*betaz*sgnbx;
  r2[5]=+sro*betay*sgnbx;
  r2[6]=0;
  r6[0]=+r2[0];
  r6[1]=+r2[1];
  r6[2]=-r2[2];
  r6[3]=-r2[3];
  r6[4]=+r2[4];
  r6[5]=+r2[5];
  r6[6]=+r2[6];
  r3[0]=ro*alphas;
  r3[1]=-csals;
  r3[2]=-cfalf*sgnbx*betay;
  r3[3]=-cfalf*sgnbx*betaz;
  r3[4]=-aaalf*sro*betay;
  r3[5]=-aaalf*sro*betaz;
  r3[6]=alphas*gamma*pr;
  r5[0]=+r3[0];
  r5[1]=-r3[1];
  r5[2]=-r3[2];
  r5[3]=-r3[3];
  r5[4]=+r3[4];
  r5[5]=+r3[5];
  r5[6]=+r3[6];
  r4[0]=1.0;
  r4[1]=0;
  r4[2]=0;
  r4[3]=0;
  r4[4]=0;
  r4[5]=0;
  r4[6]=0;
}

inline void mhd_eigen_ae(double ro, double vx, double vy, double vz, double bx, double by, double bz, double pr, double gamma,
			 double *al1, double *al2, double *al3, double *al4, double *al5, double *al6, double *al7,
			 double *l1, double *l2, double *l3, double *l4, double *l5, double *l6, double *l7, 
			 double *r1, double *r2, double *r3, double *r4, double *r5, double *r6, double *r7)
/* Calculate approximate eigenvectors for non-conservative form */
/* Consider alfven and entropy modes, fast and slow mode is ignored */
/* Outputs: */
/* al?: Dummy. Nothing returned */
/* l?: left eigenvectors (row vectors) */
/* r?: right eigenvectors (column vectors) */
{
  /* Fast mode speed */
  double iro=1.0/ro;
  double cs2=gamma*pr*iro;
  double ca2=bx*bx*iro;
  double ct2=(by*by+bz*bz)*iro;
  double cb2=cs2+ca2+ct2;
  double cf2=0.5*(cb2+sqrt(fabs(cb2*cb2-4.0*cs2*ca2)));
  double cf=sqrt(cf2);
  
  double sro=sqrt(ro);
  double isroh=0.5/sro;
  double icf2=1.0/cf2;
  double byh=0.5*by;
  double bzh=0.5*bz;

  /* Left eigenvectors */
  l1[0]=0;
  l1[1]=+1.0;
  l1[2]=0;
  l1[3]=0;
  l1[4]=0;
  l1[5]=0;
  l1[6]=0;
  l2[0]=cf2;
  l2[1]=0;
  l2[2]=0;
  l2[3]=0;
  l2[4]=-by;
  l2[5]=-bz;
  l2[6]=-1.0;
  l3[0]=0;
  l3[1]=0;
  l3[2]=0;
  l3[3]=0;
  l3[4]=+by;
  l3[5]=+bz;
  l3[6]=+1.0;
  l4[0]=0;
  l4[1]=0;
  l4[2]=+sro;
  l4[3]=0;
  l4[4]=+1.0;
  l4[5]=0;
  l4[6]=0;
  l5[0]=0;
  l5[1]=0;
  l5[2]=+sro;
  l5[3]=0;
  l5[4]=-1.0;
  l5[5]=0;
  l5[6]=0;
  l6[0]=0;
  l6[1]=0;
  l6[2]=0;
  l6[3]=+sro;
  l6[4]=0;
  l6[5]=+1.0;
  l6[6]=0;
  l7[0]=0;
  l7[1]=0;
  l7[2]=0;
  l7[3]=+sro;
  l7[4]=0;
  l7[5]=-1.0;
  l7[6]=0;

  /* Right eigenvectors */
  r1[0]=0;
  r1[1]=+1.0;
  r1[2]=0;
  r1[3]=0;
  r1[4]=0;
  r1[5]=0;
  r1[6]=0;
  r2[0]=+icf2;
  r2[1]=0;
  r2[2]=0;
  r2[3]=0;
  r2[4]=0;
  r2[5]=0;
  r2[6]=0;
  r3[0]=+icf2;
  r3[1]=0;
  r3[2]=0;
  r3[3]=0;
  r3[4]=0;
  r3[5]=0;
  r3[6]=+1.0;
  r4[0]=0;
  r4[1]=0;
  r4[2]=+isroh;
  r4[3]=0;
  r4[4]=+0.5;
  r4[5]=0;
  r4[6]=-byh;
  r5[0]=0;
  r5[1]=0;
  r5[2]=+isroh;
  r5[3]=0;
  r5[4]=-0.5;
  r5[5]=0;
  r5[6]=+byh;
  r6[0]=0;
  r6[1]=0;
  r6[2]=0;
  r6[3]=+isroh;
  r6[4]=0;
  r6[5]=+0.5;
  r6[6]=-bzh;
  r7[0]=0;
  r7[1]=0;
  r7[2]=0;
  r7[3]=+isroh;
  r7[4]=0;
  r7[5]=-0.5;
  r7[6]=+bzh;
}

inline void mhd_c_reconst(const double *ro, const double *vx, const double *vy, const double *vz,
			  const double *by, const double *bz, const double *pr,
			  double bx, double gamma, const int s,
			  double *ul, double *ur,
			  void (*func_lr)(const double*, double*, double*))
/* Reconstruction of MHD variables with characteristic decompisition (based on ATHENA) */
/* ro ~ pr: primitive variables in the stencil */
/* bx, gamma: normal magnetic field and specific heat ratio */
/* s: number of grids in the stencil */
/* ul: primitive variables at the left-face of right boundary */
/* ur: primitive variables at the right-face of left boundary */
/* func_lr: reconstruction function with degree of freedom = s */
{
  int m,n;
  double w[7][s],wl[7],wr[7];
  double al0,al1,al2,al3,al4,al5,al6; /* Not used */
  double l0[7]={1.0,0,0,0,0,0,0};
  double l1[7]={0,1.0,0,0,0,0,0};
  double l2[7]={0,0,1.0,0,0,0,0};
  double l3[7]={0,0,0,1.0,0,0,0};
  double l4[7]={0,0,0,0,1.0,0,0};
  double l5[7]={0,0,0,0,0,1.0,0};
  double l6[7]={0,0,0,0,0,0,1.0};
  double r0[7]={1.0,0,0,0,0,0,0};
  double r1[7]={0,1.0,0,0,0,0,0};
  double r2[7]={0,0,1.0,0,0,0,0};
  double r3[7]={0,0,0,1.0,0,0,0};
  double r4[7]={0,0,0,0,1.0,0,0};
  double r5[7]={0,0,0,0,0,1.0,0};
  double r6[7]={0,0,0,0,0,0,1.0};

  /* Calculate eigenvalues and eigenvectors based on ATHENA */
  mhd_eigen_st(ro[s/2],vx[s/2],vy[s/2],vz[s/2],bx,by[s/2],bz[s/2],pr[s/2],gamma,
	       &al0,&al1,&al2,&al3,&al4,&al5,&al6,
	       &l0[0],&l1[0],&l2[0],&l3[0],&l4[0],&l5[0],&l6[0],
	       &r0[0],&r1[0],&r2[0],&r3[0],&r4[0],&r5[0],&r6[0]);

  /* Characteristic variables in stencil */
#pragma simd
  for (m=0;m<s;m++){
    w[0][m]=(+l0[0]*ro[m]+l0[1]*vx[m]+l0[2]*vy[m]+l0[3]*vz[m]+l0[4]*by[m]+l0[5]*bz[m]+l0[6]*pr[m]);
    w[1][m]=(+l1[0]*ro[m]+l1[1]*vx[m]+l1[2]*vy[m]+l1[3]*vz[m]+l1[4]*by[m]+l1[5]*bz[m]+l1[6]*pr[m]);
    w[2][m]=(+l2[0]*ro[m]+l2[1]*vx[m]+l2[2]*vy[m]+l2[3]*vz[m]+l2[4]*by[m]+l2[5]*bz[m]+l2[6]*pr[m]);
    w[3][m]=(+l3[0]*ro[m]+l3[1]*vx[m]+l3[2]*vy[m]+l3[3]*vz[m]+l3[4]*by[m]+l3[5]*bz[m]+l3[6]*pr[m]);
    w[4][m]=(+l4[0]*ro[m]+l4[1]*vx[m]+l4[2]*vy[m]+l4[3]*vz[m]+l4[4]*by[m]+l4[5]*bz[m]+l4[6]*pr[m]);
    w[5][m]=(+l5[0]*ro[m]+l5[1]*vx[m]+l5[2]*vy[m]+l5[3]*vz[m]+l5[4]*by[m]+l5[5]*bz[m]+l5[6]*pr[m]);
    w[6][m]=(+l6[0]*ro[m]+l6[1]*vx[m]+l6[2]*vy[m]+l6[3]*vz[m]+l6[4]*by[m]+l6[5]*bz[m]+l6[6]*pr[m]);
  }

  /* Reconstruction */
#pragma simd
  for (n=0;n<7;n++) func_lr(&w[n][0],&wl[n],&wr[n]);
  
  /* Primitive variables (ro,vx,vy,vz,by,bz,pr) at cell faces */
#pragma simd
  for (n=0;n<7;n++)  /* Left-face at i+1/2 */
    ul[n]=r0[n]*wl[0]+r1[n]*wl[1]+r2[n]*wl[2]+r3[n]*wl[3]+r4[n]*wl[4]+r5[n]*wl[5]+r6[n]*wl[6];
#pragma simd
  for (n=0;n<7;n++)  /* Right-face at i-1/2 */
    ur[n]=r0[n]*wr[0]+r1[n]*wr[1]+r2[n]*wr[2]+r3[n]*wr[3]+r4[n]*wr[4]+r5[n]*wr[5]+r6[n]*wr[6];

  /* To exactly keep v=0 and B=0 */
  int flg[5]={1,1,1,1,1};
  for (m=0;m<s;m++){
    if (vx[m] != 0.0) flg[0]=0;
    if (vy[m] != 0.0) flg[1]=0;
    if (vz[m] != 0.0) flg[2]=0;
    if (by[m] != 0.0) flg[3]=0;
    if (bz[m] != 0.0) flg[4]=0;
  }
  if (flg[0]) ul[1]=ur[1]=0.0; /* vx */
  if (flg[1]) ul[2]=ur[2]=0.0; /* vy */
  if (flg[2]) ul[3]=ur[3]=0.0; /* vz */
  if (flg[3]) ul[4]=ur[4]=0.0; /* By */
  if (flg[4]) ul[5]=ur[5]=0.0; /* Bz */
}

inline void mhd_a_reconst(const double *ro, const double *vx, const double *vy, const double *vz,
			  const double *by, const double *bz, const double *pr,
			  double bx, double gamma, const int s,
			  double *ul, double *ur,
			  void (*func_lr)(const double*, double*, double*))
/* Reconstruction of MHD variables with alfven characteristic decompisition */
/* ro ~ pr: primitive variables in the stencil */
/* bx, gamma: normal magnetic field and specific heat ratio */
/* s: number of grids in the stencil */
/* ul: primitive variables at the left-face of right boundary */
/* ur: primitive variables at the right-face of left boundary */
/* func_lr: reconstruction function with degree of freedom = s */
{
  int m,n;
  double w[7][s],wl[7],wr[7];
  double al0,al1,al2,al3,al4,al5,al6; /* Not used */
  double l0[7]={1.0,0,0,0,0,0,0};
  double l1[7]={0,1.0,0,0,0,0,0};
  double l2[7]={0,0,1.0,0,0,0,0};
  double l3[7]={0,0,0,1.0,0,0,0};
  double l4[7]={0,0,0,0,1.0,0,0};
  double l5[7]={0,0,0,0,0,1.0,0};
  double l6[7]={0,0,0,0,0,0,1.0};
  double r0[7]={1.0,0,0,0,0,0,0};
  double r1[7]={0,1.0,0,0,0,0,0};
  double r2[7]={0,0,1.0,0,0,0,0};
  double r3[7]={0,0,0,1.0,0,0,0};
  double r4[7]={0,0,0,0,1.0,0,0};
  double r5[7]={0,0,0,0,0,1.0,0};
  double r6[7]={0,0,0,0,0,0,1.0};

  /* Calculate alfven eigenvalues and eigenvectors */
  mhd_eigen_ae(ro[s/2],vx[s/2],vy[s/2],vz[s/2],bx,by[s/2],bz[s/2],pr[s/2],gamma,
	       &al0,&al1,&al2,&al3,&al4,&al5,&al6,
	       &l0[0],&l1[0],&l2[0],&l3[0],&l4[0],&l5[0],&l6[0],
	       &r0[0],&r1[0],&r2[0],&r3[0],&r4[0],&r5[0],&r6[0]);

  /* Characteristic variables in stencil */
#pragma simd
  for (m=0;m<s;m++){
    w[0][m]=(+l0[0]*ro[m]+l0[1]*vx[m]+l0[2]*vy[m]+l0[3]*vz[m]+l0[4]*by[m]+l0[5]*bz[m]+l0[6]*pr[m]);
    w[1][m]=(+l1[0]*ro[m]+l1[1]*vx[m]+l1[2]*vy[m]+l1[3]*vz[m]+l1[4]*by[m]+l1[5]*bz[m]+l1[6]*pr[m]);
    w[2][m]=(+l2[0]*ro[m]+l2[1]*vx[m]+l2[2]*vy[m]+l2[3]*vz[m]+l2[4]*by[m]+l2[5]*bz[m]+l2[6]*pr[m]);
    w[3][m]=(+l3[0]*ro[m]+l3[1]*vx[m]+l3[2]*vy[m]+l3[3]*vz[m]+l3[4]*by[m]+l3[5]*bz[m]+l3[6]*pr[m]);
    w[4][m]=(+l4[0]*ro[m]+l4[1]*vx[m]+l4[2]*vy[m]+l4[3]*vz[m]+l4[4]*by[m]+l4[5]*bz[m]+l4[6]*pr[m]);
    w[5][m]=(+l5[0]*ro[m]+l5[1]*vx[m]+l5[2]*vy[m]+l5[3]*vz[m]+l5[4]*by[m]+l5[5]*bz[m]+l5[6]*pr[m]);
    w[6][m]=(+l6[0]*ro[m]+l6[1]*vx[m]+l6[2]*vy[m]+l6[3]*vz[m]+l6[4]*by[m]+l6[5]*bz[m]+l6[6]*pr[m]);
  }

  /* Reconstruction */
#pragma simd
  for (n=0;n<7;n++) func_lr(&w[n][0],&wl[n],&wr[n]);
  
  /* Primitive variables (ro,vx,vy,vz,by,bz,pr) at cell faces */
#pragma simd
  for (n=0;n<7;n++)  /* Left-face at i+1/2 */
    ul[n]=r0[n]*wl[0]+r1[n]*wl[1]+r2[n]*wl[2]+r3[n]*wl[3]+r4[n]*wl[4]+r5[n]*wl[5]+r6[n]*wl[6];
#pragma simd
  for (n=0;n<7;n++)  /* Right-face at i-1/2 */
    ur[n]=r0[n]*wr[0]+r1[n]*wr[1]+r2[n]*wr[2]+r3[n]*wr[3]+r4[n]*wr[4]+r5[n]*wr[5]+r6[n]*wr[6];

  /* To exactly keep v=0 and B=0 */
  int flg[5]={1,1,1,1,1};
  for (m=0;m<s;m++){
    if (vx[m] != 0.0) flg[0]=0;
    if (vy[m] != 0.0) flg[1]=0;
    if (vz[m] != 0.0) flg[2]=0;
    if (by[m] != 0.0) flg[3]=0;
    if (bz[m] != 0.0) flg[4]=0;
  }
  if (flg[0]) ul[1]=ur[1]=0.0; /* vx */
  if (flg[1]) ul[2]=ur[2]=0.0; /* vy */
  if (flg[2]) ul[3]=ur[3]=0.0; /* vz */
  if (flg[3]) ul[4]=ur[4]=0.0; /* By */
  if (flg[4]) ul[5]=ur[5]=0.0; /* Bz */
}

inline void mhd_reconst(const double *ro, const double *vx, const double *vy, const double *vz,
			const double *by, const double *bz, const double *pr,
			double bx, double gamma, const int s,
			double *ul, double *ur,
			void (*func_lr)(const double*, double*, double*))
/* Reconstruction of MHD variables  */
/* Caracteristic decomposition NOT used. Mainly used for debug */
/* ro ~ pr: primitive variables in the stencil */
/* bx, gamma: normal magnetic field and specific heat ratio */
/* s: number of grids in the stencil */
/* ul: primitive variables at the left-face of right boundary */
/* ur: primitive variables at the right-face of left boundary */
/* func_lr: reconstruction function with degree of freedom = s */
{
  int m,n;
  double w[7][s],wl[7],wr[7];

  /* Characteristic variables in stencil */
#pragma simd
  for (m=0;m<s;m++){
    w[0][m]=ro[m];
    w[1][m]=vx[m];
    w[2][m]=vy[m];
    w[3][m]=vz[m];
    w[4][m]=by[m];
    w[5][m]=bz[m];
    w[6][m]=pr[m];
  }

  /* Reconstruction */
#pragma simd
  for (n=0;n<7;n++) func_lr(&w[n][0],&wl[n],&wr[n]);
  
  /* Primitive variables (ro,vx,vy,vz,by,bz,pr) at cell faces */
  #pragma simd
  for (n=0;n<7;n++){
    ul[n]=wl[n];		/* Left-face at i+1/2 */
    ur[n]=wr[n];		/* Right-face at i-1/2 */
  }
}

#endif
