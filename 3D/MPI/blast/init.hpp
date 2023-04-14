void init_grid(int mpi_rank);
void init_plasma(int mpi_rank);

void init_grid(int mpi_rank)
// Define X,Y, and Z coordinates
{
  int m_xy=mpi_numx*mpi_numy;
  for (int i=0;i<nx;i++) x[i]=(i-XOFF+0.5+((mpi_rank%m_xy)%mpi_numx)*(XMESH/MNP_X))*dx+xmin;
  for (int j=0;j<ny;j++) y[j]=(j-YOFF+0.5+((mpi_rank%m_xy)/mpi_numx)*(YMESH/MNP_Y))*dy+ymin;
  for (int k=0;k<nz;k++) z[k]=(k-ZOFF+0.5+(mpi_rank/m_xy)*(ZMESH/MNP_Z))*dz+zmin;
}

void init_plasma(int mpi_rank)
// Set initial condition
{
  // Blast wave
  for (int k=0;k<nz;k++){
    for (int j=0;j<ny;j++){
      for (int i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;
	double rr,vx,vy,vz,pr;
	rr=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);
      
	ro[ss]=(rr <= r_0)?ro1:ro0;
	vx=0.0;
	vy=0.0;
	vz=0.0;
	pr=(rr <= r_0)?pr1:pr0;
      
	bx[ss]=b_0*sin(bthe*dtor)*cos(bphi*dtor);
	by[ss]=b_0*sin(bthe*dtor)*sin(bphi*dtor);
	bz[ss]=b_0*cos(bthe*dtor);

	mx[ss]=ro[ss]*vx;
	my[ss]=ro[ss]*vy;
	mz[ss]=ro[ss]*vz;
	en[ss]=pr/(gam-1)+0.5*(ro[ss]*(vx*vx+vy*vy+vz*vz)+(bx[ss]*bx[ss]+by[ss]*by[ss]+bz[ss]*bz[ss]));
      }
    }
  }
}
