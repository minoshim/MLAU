void mhd_fd_ct_2d(double *ro, double *mx, double *my, double *mz,
		  double *en, double *bx, double *by, double *bz,
		  double dt, double dx, double dy,
		  int nx, int ny, int xoff, int yoff, double gamma);
void hall_fd_ct_2d(double *ro, double *mx, double *my, double *mz,
		   double *en, double *bx, double *by, double *bz,
		   double dt, double dx, double dy,
		   int nx, int ny, int xoff, int yoff, double eta_h);
void bc2d(double *f, int nx, int ny, int xoff, int yoff,
	  int stx, int dnx, int sty, int dny);

inline double fmode(double rho, double p, double bsqr, double gamma)
{
  /* Fast mode velocity in perpendicular direction */
  /* Enough for determing maximum characteristic velocity */
  return( sqrt((0.5*gamma*p+bsqr)/(rho)) );
}
