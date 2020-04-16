void init_grid(int mpi_rank)
// Define X and Y coordinates
{
  for (int i=0;i<nx;i++) x[i]=(i-XOFF+0.5+(mpi_rank%mpi_numx)*(XMESH/MNP_X))*dx;
  for (int j=0;j<ny;j++) y[j]=(j-YOFF+0.5+(mpi_rank/mpi_numx)*(YMESH/MNP_Y))*dy;
}
