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

void mpi_ybc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D Y Dirichlet or Neumann BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1) or Neumann (+1). if dn==0, nothing to do */
{
  if (dn != 0){
    int i,j;
    /* Left */
    if (mpi_rank/mpi_numx == 0){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++) f[nx*j+i]=dn*f[nx*((2*yoff-1+st)-j)+i];
      }
    }
    /* Right */
    if (mpi_rank/mpi_numx == (mpi_numy-1)){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++) f[nx*(ny-1-j)+i]=dn*f[nx*((ny-2*yoff+st)+j)+i];
      }
    }
  }
}
