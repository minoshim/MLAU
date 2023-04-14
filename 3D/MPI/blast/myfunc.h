void mpi_sdrv3d(double *f[], int nn, int nx, int ny, int nz, int xoff, int yoff, int zoff,
		int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);
void mpi_xbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);
void mpi_ybc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);
void mpi_zbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);
void mhd_fd_ct_3d(double *ro, double *mx, double *my, double *mz,
		  double *en, double *bx, double *by, double *bz,
		  double dt, double dx, double dy, double dz,
		  int nx, int ny, int nz, int xoff, int yoff, int zoff, double gamma,
		  int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);

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

