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

void mpi_xbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);
void mpi_ybc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);
void mpi_zbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);

void mpi_xbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D X Dirichlet or Neumann BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1) or Neumann (+1). if dn==0, nothing to do */
{
  if (dn != 0){
    unsigned long i,j,k;
    int m_xy=mpi_numx*mpi_numy;
    /* Left */
    if (((mpi_rank%m_xy)%mpi_numx) == 0){
      for (k=0;k<nz;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<xoff;i++) f[nx*(ny*k+j)+i]=dn*f[nx*(ny*k+j)+(2*xoff-1+st)-i];
	}
      }
    }
    /* Right */
    if (((mpi_rank%m_xy)%mpi_numx) == (mpi_numx-1)){
      for (k=0;k<nz;k++){
        for (j=0;j<ny;j++){
	  for (i=0;i<xoff;i++) f[nx*(ny*k+j)+(nx-1-i)]=dn*f[nx*(ny*k+j)+(nx-2*xoff+st)+i];
        }
      }
    }
  }
}

void mpi_ybc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D Y Dirichlet or Neumann BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1) or Neumann (+1). if dn==0, nothing to do */
{
  if (dn != 0){
    unsigned long i,j,k;
    int m_xy=mpi_numx*mpi_numy;
    /* Left */
    if (((mpi_rank%m_xy)/mpi_numx) == 0){
      for (k=0;k<nz;k++){
	for (j=0;j<yoff;j++){
	  for (i=0;i<nx;i++) f[nx*(ny*k+j)+i]=dn*f[nx*(ny*k+(2*yoff-1+st)-j)+i];
	}
      }
    }
    /* Right */
    if (((mpi_rank%m_xy)/mpi_numx) == (mpi_numy-1)){
      for (k=0;k<nz;k++){
	for (j=0;j<yoff;j++){
	  for (i=0;i<nx;i++) f[nx*(ny*k+(ny-1-j))+i]=dn*f[nx*(ny*k+(ny-2*yoff+st)+j)+i];
	}
      }
    }
  }
}

void mpi_zbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D Z Dirichlet or Neumann BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1) or Neumann (+1). if dn==0, nothing to do */
{
  if (dn != 0){
    unsigned long i,j,k;
    int m_xy=mpi_numx*mpi_numy;
    /* Left */
    if ((mpi_rank/m_xy) == 0){
      for (k=0;k<zoff;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<nx;i++) f[nx*(ny*k+j)+i]=dn*f[nx*(ny*((2*zoff-1+st)-k)+j)+i];
	}
      }
    }
    /* Right */
    if ((mpi_rank/m_xy) == (mpi_numz-1)){
      for (k=0;k<zoff;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<nx;i++) f[nx*(ny*(nz-1-k)+j)+i]=dn*f[nx*(ny*((nz-2*zoff+st)+k)+j)+i];
	}
      }
    }
  }
}
