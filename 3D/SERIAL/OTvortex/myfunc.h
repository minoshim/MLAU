void bc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff,
	  int stx, int dnx, int sty, int dny, int stz, int dnz);

void mhd_fd_ct_3d(double *ro, double *mx, double *my, double *mz,
		  double *en, double *bx, double *by, double *bz,
		  double dt, double dx, double dy, double dz,
		  int nx, int ny, int nz, int xoff, int yoff, int zoff, double gamma);

inline double min3(double a, double b, double c)
{
  double ans=a;
  if (b < ans) ans=b;
  if (c < ans) ans=c;
  return(ans);
}
inline double fmode(double rho, double p, double bsqr, double gamma)
{
  /* Fast mode velocity in perpendicular direction */
  /* Enough for determing maximum characteristic velocity */
  return( sqrt((0.5*gamma*p+bsqr)/(rho)) );
}

