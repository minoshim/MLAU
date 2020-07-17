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

#ifndef _INTERP_H_
#define _INTERP_H_

#include <math.h>
#include "funcs.h"

/* MUSCL routines */
inline double muscl_mm_cal_f(double fu, double f0, double fd)
{
  /* Return f between f0 and fd with MinMod limiter */
  return(f0-0.5*minmod(fu-f0,f0-fd));
}
inline void muscl_mm_cal_flr(double f1, double f2, double f3,
			     double *fl, double *fr)
{
  /* Return f2.5_l and f1.5_r with MinMod limiter */
  double df=0.5*minmod(f2-f1,f3-f2);
  (*fl)=f2+df;
  (*fr)=f2-df;
}

inline double muscl_mc_cal_f(double fu, double f0, double fd)
{
  /* Return f between f0 and fd with MC limiter */
  double alpha=2.0;
  return(f0-0.5*minmod(alpha*minmod(fu-f0,f0-fd),0.5*(fu-fd)));
}
inline void muscl_mc_cal_flr(double f1, double f2, double f3,
			     double *fl, double *fr)
{
  /* Return f2.5_l and f1.5_r with MC limiter */
  double alpha=2.0;
  double df=0.5*minmod(alpha*minmod(f2-f1,f3-f2),0.5*(f3-f1));
  (*fl)=f2+df;
  (*fr)=f2-df;
}

inline double muscl_kr_cal_f(double fu, double f0, double fd)
{
  /* Return f between f0 and fd with Koren limiter */
  double df1=2.0*(fu-f0);
  double df2=2.0*(f0-fd);
  return(f0-0.5*minmod3(df1,df2,(0.5*df1+df2)/3.0));
}
inline void muscl_kr_cal_flr(double f1, double f2, double f3,
			     double *fl, double *fr)
{
  /* Return f2.5_l and f1.5_r with Koren limiter */
  double df1=2.0*(f2-f1);
  double df2=2.0*(f3-f2);
  double dfl=0.5*minmod3(df1,df2,(0.5*df1+df2)/3.0);
  double dfr=0.5*minmod3(df1,df2,(df1+0.5*df2)/3.0);
  (*fl)=f2+dfl;
  (*fr)=f2-dfr;
}

/* WCNS routines */
inline double wcns3_cal_f(double fup, double f, double fdn)
/* 3rd order WCNS interpolation between f and fdn */
/* Use ESWENO3 modification by Yamaleev+09 */
/* Minoshima+19, ApJS,242,14 */
{
  double is[2],w[2],c[2];
  double eps=1e-40,c1,tau3;

  /* Smoothness indicator */
  c1=fup-f;
  is[0]=c1*c1;
  c1=fdn-f;
  is[1]=c1*c1;

  /* Yamaleev+ 09 */
  tau3=fup-2.0*f+fdn;
  tau3=tau3*tau3;

  /* Linear weight (3rd order) */
  c[0]=0.25;
  c[1]=0.75;

  /* Nonlinear weight */
  w[0]=c[0]*(1.0+tau3/(is[0]+eps));
  w[1]=c[1]*(1.0+tau3/(is[1]+eps));
  w[0]/=(w[0]+w[1]);
  w[1]=1.0-w[0];

  return( 0.5*(w[0]*(-fup+3.0*f)+w[1]*(f+fdn)) );
}

inline void wcns3_cal_flr(double f0, double f1, double f2, double *fl, double *fr)
/* 3rd order WCNS interpolation */
/* All f are point values */
/* fl = Left-side value at 1.5 */
/* fr = Right-side value at 0.5 */
/* Minoshima+19, ApJS,242,14 */
{
  double is[2],w[2],c[2];
  double eps=1e-40,c1,tau3;

  /* Smoothness indicator */
  c1=f0-f1;
  is[0]=c1*c1;
  c1=f2-f1;
  is[1]=c1*c1;

  /* Yamaleev+ 09 */
  tau3=f2-2.0*f1+f0;
  tau3=tau3*tau3;

  is[0]=(1.0+tau3/(is[0]+eps));
  is[1]=(1.0+tau3/(is[1]+eps));

  /* Linear weight (3rd order) */
  c[0]=0.25;
  c[1]=0.75;

  /* Left-side value at 1.5 */
  w[0]=c[0]*is[0];
  w[1]=c[1]*is[1];
  w[0]/=(w[0]+w[1]);
  w[1]=1.0-w[0];
  *fl=0.5*(w[0]*(-f0+3.0*f1)+w[1]*(f1+f2));

  /* Right-side value at 0.5 */
  w[0]=c[1]*is[0];
  w[1]=c[0]*is[1];
  w[0]/=(w[0]+w[1]);
  w[1]=1.0-w[0];
  *fr=0.5*(w[1]*(-f2+3.0*f1)+w[0]*(f1+f0));
}

inline double wcns4_cal_f(double fuu, double fup, double f, double fdn, double fdd)
/* 4th order WCNS interpolation between f and fdn */
/* All f are point values */
/* Minoshima+19, ApJS,242,14 */
{
  double is[3],w[3],c[3];
  double eps=1e-40,c1,c2,denom;

  /* Smoothness indicator */
  c1=fuu-4.0*fup+3.0*f;
  c2=fuu-2.0*fup+f;
  is[0]=(13./12.)*c2*c2+0.25*c1*c1;
  c1=fup-fdn;
  c2=fup-2.0*f+fdn;
  is[1]=(13./12.)*c2*c2+0.25*c1*c1;
  c1=3.0*f-4.0*fdn+fdd;
  c2=f-2.0*fdn+fdd;
  is[2]=(13./12.)*c2*c2+0.25*c1*c1;

  /* Linear weight, achieving 5th order of dfdx=(27*(fi+1/2-fi-1/2)-(fi+3/2-fi-3/2))/24 */
  c[0]=0.075;
  c[1]=0.650;
  c[2]=0.275;

  /* Nonlinear weight */
  w[0]=c[0]/(is[0]*is[0]+eps);
  w[1]=c[1]/(is[1]*is[1]+eps);
  w[2]=c[2]/(is[2]*is[2]+eps);
  denom=1.0/(w[0]+w[1]+w[2]);
  w[0]*=denom;
  w[1]*=denom;
  w[2]=1.0-w[0]-w[1];

  return(0.125*(+w[0]*(3.0*fuu-10.0*fup+15.0*f)
		+w[1]*(-fup+6.0*f+3.0*fdn)
		+w[2]*(+3.0*f+6.0*fdn-fdd)));
}

inline void wcns4_cal_flr(double f0, double f1, double f2, double f3, double f4, double *fl, double *fr)
/* 4th order WCNS interpolatoin */
/* All f are point values */
/* fl = Left-side value at 2.5 */
/* fr = Right-side value at 1.5 */
/* Minoshima+19, ApJS,242,14 */
{
  double is[3],w[3],c[3];
  double eps=1e-40,c1,c2,denom;

  /* Smoothness indicator */
  c1=f0-4.0*f1+3.0*f2;
  c2=f0-2.0*f1+f2;
  is[0]=(13./12.)*c2*c2+0.25*c1*c1;
  c1=f1-f3;
  c2=f1-2.0*f2+f3;
  is[1]=(13./12.)*c2*c2+0.25*c1*c1;
  c1=f4-4.0*f3+3.0*f2;
  c2=f4-2.0*f3+f2;
  is[2]=(13./12.)*c2*c2+0.25*c1*c1;
  is[0]=1.0/(is[0]*is[0]+eps);
  is[1]=1.0/(is[1]*is[1]+eps);
  is[2]=1.0/(is[2]*is[2]+eps);

  /* Linear weight, achieving 5th order of dfdx=(27*(fi+1/2-fi-1/2)-(fi+3/2-fi-3/2))/24 */
  c[0]=0.075;
  c[1]=0.650;
  c[2]=0.275;

  /* Left-side value at 2.5 */
  w[0]=c[0]*is[0];
  w[1]=c[1]*is[1];
  w[2]=c[2]*is[2];
  denom=1.0/(w[0]+w[1]+w[2]);
  w[0]*=denom;
  w[1]*=denom;
  w[2]=1.0-w[0]-w[1];
  *fl=0.125*(+w[0]*(3.0*f0-10.0*f1+15.0*f2)
	     +w[1]*(-f1+6.0*f2+3.0*f3)
	     +w[2]*(+3.0*f2+6.0*f3-f4));

  /* Right-side value at 1.5 */
  w[0]=c[2]*is[0];
  w[1]=c[1]*is[1];
  w[2]=c[0]*is[2];
  denom=1.0/(w[0]+w[1]+w[2]);
  w[0]*=denom;
  w[1]*=denom;
  w[2]=1.0-w[0]-w[1];
  *fr=0.125*(+w[2]*(3.0*f4-10.0*f3+15.0*f2)
	     +w[1]*(-f3+6.0*f2+3.0*f1)
	     +w[0]*(+3.0*f2+6.0*f1-f0));
}

#endif
