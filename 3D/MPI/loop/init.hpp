#define ODR (2)			// Set 4 if high order (>2) MHD scheme is used.

void init_grid(int mpi_rank);
void init_plasma(int mpi_rank);

inline double cal_r(double x, double y, double z, double cost, double sint)
{
  double x1=+x*cost+z*sint;
  double x2=y;
  double x3=-x*sint+z*cost;
  // Satisfy periodic boundary condition (see gardiner+08)
  while(x1 > +0.5*cost){
    x1-=cost;
  }
  while(x1 < -0.5*cost){
    x1+=cost;
  }
  while(x2 > +0.5){
    x2-=1.0;
  }
  while(x2 < -0.5){
    x2+=1.0;
  }
  return(sqrt(x1*x1+x2*x2));
}

inline double cal_a3(double r, double rad, double a0)
{
  double ans=0;
  if (r <= rad) ans=a0*(rad-r);
  return(ans);
}

inline double cal_da(double *a)
{
  double ans=(a[2]-a[1]);	// 2nd order discretization
#if (ODR > 2)
  ans=(27.0*(a[2]-a[1])-(a[3]-a[0]))/24.0; // 4th order discretization
#endif
  return(ans);
}

void init_grid(int mpi_rank)
// Define X,Y, and Z coordinates
{
  int m_xy=mpi_numx*mpi_numy;
  for (int i=0;i<nx;i++) x[i]=(i-xoff+0.5+((mpi_rank%m_xy)%mpi_numx)*(XMESH/MNP_X))*dx+xmin;
  for (int j=0;j<ny;j++) y[j]=(j-yoff+0.5+((mpi_rank%m_xy)/mpi_numx)*(YMESH/MNP_Y))*dy+ymin;
  for (int k=0;k<nz;k++) z[k]=(k-zoff+0.5+(mpi_rank/m_xy)*(ZMESH/MNP_Z))*dz+zmin;
}

void init_plasma(int mpi_rank)
// Set initial condition
{
  double cost=cos(tilt);
  double sint=sin(tilt);
  
  // 3D loop advection based on Lee+13
  int i,j,k,ss;
  for (k=0;k<nz;k++){
    for (j=0;j<ny;j++){
      for (i=0;i<nx;i++){
	double vx,vy,vz;	
	ss=nx*(ny*k+j)+i;

	ro[ss]=ro0;
	vx=vx0;
	vy=vy0;
	vz=vz0;
	mx[ss]=ro[ss]*vx;
	my[ss]=ro[ss]*vy;
	mz[ss]=ro[ss]*vz;

	double rr[4],a3[4];
	// CT Bx
	rr[0]=cal_r(x[i]-0.5*dx,y[j]-1.5*dy,z[k],cost,sint);
	rr[1]=cal_r(x[i]-0.5*dx,y[j]-0.5*dy,z[k],cost,sint);
	rr[2]=cal_r(x[i]-0.5*dx,y[j]+0.5*dy,z[k],cost,sint);
	rr[3]=cal_r(x[i]-0.5*dx,y[j]+1.5*dy,z[k],cost,sint);
	a3[0]=cal_a3(rr[0],rad,a0);
	a3[1]=cal_a3(rr[1],rad,a0);
	a3[2]=cal_a3(rr[2],rad,a0);
	a3[3]=cal_a3(rr[3],rad,a0);
	bx[ss]=+cost*idy*cal_da(a3);
	// CT By
	rr[0]=cal_r(x[i]-1.5*dx,y[j]-0.5*dy,z[k],cost,sint);
	rr[1]=cal_r(x[i]-0.5*dx,y[j]-0.5*dy,z[k],cost,sint);
	rr[2]=cal_r(x[i]+0.5*dx,y[j]-0.5*dy,z[k],cost,sint);
	rr[3]=cal_r(x[i]+1.5*dx,y[j]-0.5*dy,z[k],cost,sint);
	a3[0]=cal_a3(rr[0],rad,a0);
	a3[1]=cal_a3(rr[1],rad,a0);
	a3[2]=cal_a3(rr[2],rad,a0);
	a3[3]=cal_a3(rr[3],rad,a0);
	by[ss]=-cost*idx*cal_da(a3);
	rr[0]=cal_r(x[i],y[j]-0.5*dy,z[k]-1.5*dz,cost,sint);
	rr[1]=cal_r(x[i],y[j]-0.5*dy,z[k]-0.5*dz,cost,sint);
	rr[2]=cal_r(x[i],y[j]-0.5*dy,z[k]+0.5*dz,cost,sint);
	rr[3]=cal_r(x[i],y[j]-0.5*dy,z[k]+1.5*dz,cost,sint);
	a3[0]=cal_a3(rr[0],rad,a0);
	a3[1]=cal_a3(rr[1],rad,a0);
	a3[2]=cal_a3(rr[2],rad,a0);
	a3[3]=cal_a3(rr[3],rad,a0);
	by[ss]-=sint*idz*cal_da(a3);
	// CT Bz
	rr[0]=cal_r(x[i],y[j]-1.5*dy,z[k]-0.5*dz,cost,sint);
	rr[1]=cal_r(x[i],y[j]-0.5*dy,z[k]-0.5*dz,cost,sint);
	rr[2]=cal_r(x[i],y[j]+0.5*dy,z[k]-0.5*dz,cost,sint);
	rr[3]=cal_r(x[i],y[j]+1.5*dy,z[k]-0.5*dz,cost,sint);
	a3[0]=cal_a3(rr[0],rad,a0);
	a3[1]=cal_a3(rr[1],rad,a0);
	a3[2]=cal_a3(rr[2],rad,a0);
	a3[3]=cal_a3(rr[3],rad,a0);
	bz[ss]=+sint*idy*cal_da(a3);
      }
    }
  }

  // Energy
  for (k=1;k<nz-1;k++){
    for (j=1;j<ny-1;j++){
      for (i=1;i<nx-1;i++){
	double pr,cx,cy,cz;
	ss=nx*(ny*k+j)+i;
	pr=pr0;
	cx=0.5*(bx[ss]+bx[nx*(ny*k+j)+(i+1)]);
	cy=0.5*(by[ss]+by[nx*(ny*k+(j+1))+i]);
	cz=0.5*(bz[ss]+bz[nx*(ny*(k+1)+j)+i]);
	en[ss]=pr/(gam-1)+0.5*((mx[ss]*mx[ss]+my[ss]*my[ss]+mz[ss]*mz[ss])/ro[ss]+(cx*cx+cy*cy+cz*cz));
      }
    }
  }
}
