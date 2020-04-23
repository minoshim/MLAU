void init_grid(int mpi_rank)
// Define X,Y, and Z coordinates
{
  int m_xy=mpi_numx*mpi_numy;
  for (int i=0;i<nx;i++) x[i]=(i-XOFF+((mpi_rank%m_xy)%mpi_numx)*(XMESH/MNP_X))*dx;
  for (int j=0;j<ny;j++) y[j]=(j-YOFF+((mpi_rank%m_xy)/mpi_numx)*(YMESH/MNP_Y))*dy;
  for (int k=0;k<nz;k++) z[k]=(k-ZOFF+(mpi_rank/m_xy)*(ZMESH/MNP_Z))*dz;
}
