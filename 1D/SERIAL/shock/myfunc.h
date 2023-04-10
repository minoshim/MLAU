#ifndef _MYFUNC_H_
#define _MYFUNC_H_

#include "funcs.h"

inline double sqrwave(double x, double x0, double x1, double dx)
{
  /* Return square wave distribution */
  return(0.5*(tanh((x-x0)/dx)-tanh((x-x1)/dx)));
}
inline double rand_noise(const double *params, unsigned seed)
{
  /* Return uniform random distribution */
  /* Give seed */
  static int r_flag=0;
  if (r_flag == 0){
    srandom(seed);
    r_flag=1;
  }
  return(params[0]+params[1]*((double)random()/RAND_MAX-0.5)*2.0);
}

void openbc1d_mhd(double *ro, double *mx, double *my, double *mz,
		  double *en, double *by, double *bz,
		  double bx,
		  int nx, int xoff, double gamma,
		  int istt, int iend);
void mhd_fd4c_1d(double *ro, double *mx, double *my, double *mz,
		 double *en, double *by, double *bz,
		 double bx, double dt, double dx,
		 int nx, int xoff, double gamma);

#endif
