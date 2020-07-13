void mpi_sdrv2d(double *f[], int nn, int nx, int ny, int xoff, int yoff, 
		int mpi_rank, int mpi_numx, int mpi_numy);
void mpi_xbc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy);
void mpi_ybc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy);
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
  
inline double harris_field(double x, const double *params)
{
  /* Harris magnetic field */
  double x0=params[0],width=params[1]+1e-15;
  return( tanh((x-x0)/width) );
}
inline double harris_density(double x, const double *params)
{
  /* Harris distribution of the density */
  double x0=params[0],width=params[1]+1e-15;
  double cosh1=cosh((x-x0)/width);
  return( 1.0/(cosh1*cosh1) );
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
inline double fmode(double rho, double p, double bsqr, double gamma)
{
  /* Fast mode velocity in perpendicular direction */
  /* Enough for determing maximum characteristic velocity */
  return( sqrt((0.5*gamma*p+bsqr)/(rho)) );
}
