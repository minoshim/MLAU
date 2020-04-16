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
#include "mhd_eigen.h"

void openbc1d_mhd(double *ro, double *mx, double *my, double *mz,
		  double *en, double *by, double *bz,
		  double bx,
		  int nx, int xoff, double gamma,
		  int istt, int iend)
/* Open Boundary condition for MHD */
/* Use characteristic decomposition */

/* istt and iend are indices to set boundary condition */
/* [istt,iend]=[0,1] => set left and right boundaries */
/* [istt,iend]=[0,0] => set left boundary only */
/* [istt,iend]=[1,1] => set right boundary only */
{
  int i;
  const int bpos[2]={xoff,nx-xoff-1}; /* Edge of simulation domain */
  const int dirc[2]={1,-1};	      /* 1 for Left, -1 for Right */
  const double bx2=bx*bx;
  
  for (i=istt;i<=iend;i++){		/* Left @ i=0, Right @ i=1 */
    int sgn,si,so,ss,j;
    double iro,v2,b2;
    double rob[2],vxb[2],vyb[2],vzb[2],bxb[2],byb[2],bzb[2],prb[2];
    double w[7][2];
    double al0,al1,al2,al3,al4,al5,al6;
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

    sgn=dirc[i];    		/* 1 for Left, -1 for Right */
    si=bpos[i];			/* Position @ inside the boundary */
    so=bpos[i]-sgn;		/* Position @ outside the boundary */

    /* Primitive variables inside the boundary */
    {
      ss=si;
      j=0;
      rob[j]=ro[ss];
      iro=1.0/rob[j];
      vxb[j]=mx[ss]*iro;			/* vx */
      vyb[j]=my[ss]*iro;			/* vy */
      vzb[j]=mz[ss]*iro;			/* vz */
      bxb[j]=bx;
      byb[j]=by[ss];
      bzb[j]=bz[ss];
      v2=vxb[j]*vxb[j]+vyb[j]*vyb[j]+vzb[j]*vzb[j];
      b2=bxb[j]*bxb[j]+byb[j]*byb[j]+bzb[j]*bzb[j];
      prb[j]=(gamma-1)*(en[ss]-0.5*(rob[j]*v2+b2)); /* pr */
    }
    /* Primitive variables outside the boundary */
    {
      ss=so;
      j=1;
      rob[j]=ro[ss];
      iro=1.0/rob[j];
      vxb[j]=mx[ss]*iro;			/* vx */
      vyb[j]=my[ss]*iro;			/* vy */
      vzb[j]=mz[ss]*iro;			/* vz */
      bxb[j]=bx;
      byb[j]=by[ss];
      bzb[j]=bz[ss];
      v2=vxb[j]*vxb[j]+vyb[j]*vyb[j]+vzb[j]*vzb[j];
      b2=bxb[j]*bxb[j]+byb[j]*byb[j]+bzb[j]*bzb[j];
      prb[j]=(gamma-1)*(en[ss]-0.5*(rob[j]*v2+b2)); /* pr */
    }
    /* Convert to characteristic variables */
    {
      /* Eigenvectors @ the boundary. Use arithmetic average for simplcility */
      mhd_eigen_st(0.5*(rob[0]+rob[1]),0.5*(vxb[0]+vxb[1]),0.5*(vyb[0]+vyb[1]),0.5*(vzb[0]+vzb[1]),
		   0.5*(bxb[0]+bxb[1]),0.5*(byb[0]+byb[1]),0.5*(bzb[0]+bzb[1]),0.5*(prb[0]+prb[1]),gamma,
		   &al0,&al1,&al2,&al3,&al4,&al5,&al6,
		   &l0[0],&l1[0],&l2[0],&l3[0],&l4[0],&l5[0],&l6[0],
		   &r0[0],&r1[0],&r2[0],&r3[0],&r4[0],&r5[0],&r6[0]);
      /* Characteristic variables inside the boundary */
      j=0;
      w[0][j]=(+l0[0]*rob[j]+l0[1]*vxb[j]+l0[2]*vyb[j]+l0[3]*vzb[j]+l0[4]*byb[j]+l0[5]*bzb[j]+l0[6]*prb[j]);
      w[1][j]=(+l1[0]*rob[j]+l1[1]*vxb[j]+l1[2]*vyb[j]+l1[3]*vzb[j]+l1[4]*byb[j]+l1[5]*bzb[j]+l1[6]*prb[j]);
      w[2][j]=(+l2[0]*rob[j]+l2[1]*vxb[j]+l2[2]*vyb[j]+l2[3]*vzb[j]+l2[4]*byb[j]+l2[5]*bzb[j]+l2[6]*prb[j]);
      w[3][j]=(+l3[0]*rob[j]+l3[1]*vxb[j]+l3[2]*vyb[j]+l3[3]*vzb[j]+l3[4]*byb[j]+l3[5]*bzb[j]+l3[6]*prb[j]);
      w[4][j]=(+l4[0]*rob[j]+l4[1]*vxb[j]+l4[2]*vyb[j]+l4[3]*vzb[j]+l4[4]*byb[j]+l4[5]*bzb[j]+l4[6]*prb[j]);
      w[5][j]=(+l5[0]*rob[j]+l5[1]*vxb[j]+l5[2]*vyb[j]+l5[3]*vzb[j]+l5[4]*byb[j]+l5[5]*bzb[j]+l5[6]*prb[j]);
      w[6][j]=(+l6[0]*rob[j]+l6[1]*vxb[j]+l6[2]*vyb[j]+l6[3]*vzb[j]+l6[4]*byb[j]+l6[5]*bzb[j]+l6[6]*prb[j]);
      /* Characteristic variables outside the boundary */
      j=1;
      w[0][j]=(+l0[0]*rob[j]+l0[1]*vxb[j]+l0[2]*vyb[j]+l0[3]*vzb[j]+l0[4]*byb[j]+l0[5]*bzb[j]+l0[6]*prb[j]);
      w[1][j]=(+l1[0]*rob[j]+l1[1]*vxb[j]+l1[2]*vyb[j]+l1[3]*vzb[j]+l1[4]*byb[j]+l1[5]*bzb[j]+l1[6]*prb[j]);
      w[2][j]=(+l2[0]*rob[j]+l2[1]*vxb[j]+l2[2]*vyb[j]+l2[3]*vzb[j]+l2[4]*byb[j]+l2[5]*bzb[j]+l2[6]*prb[j]);
      w[3][j]=(+l3[0]*rob[j]+l3[1]*vxb[j]+l3[2]*vyb[j]+l3[3]*vzb[j]+l3[4]*byb[j]+l3[5]*bzb[j]+l3[6]*prb[j]);
      w[4][j]=(+l4[0]*rob[j]+l4[1]*vxb[j]+l4[2]*vyb[j]+l4[3]*vzb[j]+l4[4]*byb[j]+l4[5]*bzb[j]+l4[6]*prb[j]);
      w[5][j]=(+l5[0]*rob[j]+l5[1]*vxb[j]+l5[2]*vyb[j]+l5[3]*vzb[j]+l5[4]*byb[j]+l5[5]*bzb[j]+l5[6]*prb[j]);
      w[6][j]=(+l6[0]*rob[j]+l6[1]*vxb[j]+l6[2]*vyb[j]+l6[3]*vzb[j]+l6[4]*byb[j]+l6[5]*bzb[j]+l6[6]*prb[j]);
    }
    /* Set boundary condition */
    {
      ss=so;
      /* Boundary condition for characteristic variables, considering the sign of eigenvalues */
      /* If variable is incoming, keep initial data, else copy inner data  */
      w[0][1]=(sgn*al0 > 0)?w[0][1]:w[0][0];
      w[1][1]=(sgn*al1 > 0)?w[1][1]:w[1][0];
      w[2][1]=(sgn*al2 > 0)?w[2][1]:w[2][0];
      w[3][1]=(sgn*al3 > 0)?w[3][1]:w[3][0];
      w[4][1]=(sgn*al4 > 0)?w[4][1]:w[4][0];
      w[5][1]=(sgn*al5 > 0)?w[5][1]:w[5][0];
      w[6][1]=(sgn*al6 > 0)?w[6][1]:w[6][0];
      /* Primitive variables */
      ro[ss]=r0[0]*w[0][1]+r1[0]*w[1][1]+r2[0]*w[2][1]+r3[0]*w[3][1]+r4[0]*w[4][1]+r5[0]*w[5][1]+r6[0]*w[6][1];
      mx[ss]=r0[1]*w[0][1]+r1[1]*w[1][1]+r2[1]*w[2][1]+r3[1]*w[3][1]+r4[1]*w[4][1]+r5[1]*w[5][1]+r6[1]*w[6][1]; /* vx */
      my[ss]=r0[2]*w[0][1]+r1[2]*w[1][1]+r2[2]*w[2][1]+r3[2]*w[3][1]+r4[2]*w[4][1]+r5[2]*w[5][1]+r6[2]*w[6][1]; /* vy */
      mz[ss]=r0[3]*w[0][1]+r1[3]*w[1][1]+r2[3]*w[2][1]+r3[3]*w[3][1]+r4[3]*w[4][1]+r5[3]*w[5][1]+r6[3]*w[6][1]; /* vz */
      by[ss]=r0[4]*w[0][1]+r1[4]*w[1][1]+r2[4]*w[2][1]+r3[4]*w[3][1]+r4[4]*w[4][1]+r5[4]*w[5][1]+r6[4]*w[6][1];
      bz[ss]=r0[5]*w[0][1]+r1[5]*w[1][1]+r2[5]*w[2][1]+r3[5]*w[3][1]+r4[5]*w[4][1]+r5[5]*w[5][1]+r6[5]*w[6][1];
      en[ss]=r0[6]*w[0][1]+r1[6]*w[1][1]+r2[6]*w[2][1]+r3[6]*w[3][1]+r4[6]*w[4][1]+r5[6]*w[5][1]+r6[6]*w[6][1]; /* pr */
      /* Convert to conservative variables */
      v2=mx[ss]*mx[ss]+my[ss]*my[ss]+mz[ss]*mz[ss];
      b2=bx2+by[ss]*by[ss]+bz[ss]*bz[ss];
      en[ss]=en[ss]/(gamma-1.0)+0.5*(ro[ss]*v2+b2);
      mx[ss]*=ro[ss];
      my[ss]*=ro[ss];
      mz[ss]*=ro[ss];
    }

    /* Remaining copy */
    for (j=0;j<xoff-1;j++){
      ss=i*(nx-1)+sgn*j;
      ro[ss]=ro[so];
      mx[ss]=mx[so];
      my[ss]=my[so];
      mz[ss]=mz[so];
      by[ss]=by[so];
      bz[ss]=bz[so];
      en[ss]=en[so];
    }
  }
  
}
