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

#include <stdlib.h>
#include "mpi.h"

#define PRDC_X (1)		/* Set 1 for periodic in X */
#define PRDC_Y (1)		/* Set 1 for periodic in Y */

void mpi_sdrv2d(double *f[], int nn, int nx, int ny, int xoff, int yoff, 
		int mpi_rank, int mpi_numx, int mpi_numy);
void mpi_xbc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy);
void mpi_ybc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy);

void mpi_sdrv2d(double *f[], int nn, int nx, int ny, int xoff, int yoff, 
		int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D MPI SendRecv */
{
  int i,j,n;
  int mpi_tag=0;
  int rankl,rankh;
  int ntot,ntot2;
  MPI_Status r_stat;
  double *fold,*fcpy;

  /* XBC */
#if (PRDC_X)
  rankl=((mpi_rank % mpi_numx) == 0)?(mpi_rank+(mpi_numx-1)):(mpi_rank-1);
  rankh=((mpi_rank % mpi_numx) == (mpi_numx-1))?(mpi_rank-(mpi_numx-1)):(mpi_rank+1);
#else
  rankl=((mpi_rank % mpi_numx) == 0)?(MPI_PROC_NULL):(mpi_rank-1);
  rankh=((mpi_rank % mpi_numx) == (mpi_numx-1))?(MPI_PROC_NULL):(mpi_rank+1);
#endif
  if (mpi_numx != 1){
    ntot=nn*ny*xoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (i=0;i<xoff;i++){
      for (j=0;j<ny;j++){	/* Transpose */
	for (n=0;n<nn;n++){
	  fold[nn*(ny*i+j)+n]=f[n][nx*j+(xoff+i)];
	  fold[nn*(ny*(2*xoff-1-i)+j)+n]=f[n][nx*j+(nx-xoff-1-i)];
	  fcpy[nn*(ny*i+j)+n]=f[n][nx*j+i];
	  fcpy[nn*(ny*(2*xoff-1-i)+j)+n]=f[n][nx*j+(nx-1-i)];
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (j=0;j<ny;j++){
      for (i=0;i<xoff;i++){
	for (n=0;n<nn;n++){
	  f[n][nx*j+i]=fcpy[nn*(ny*i+j)+n];
	  f[n][nx*j+(nx-1-i)]=fcpy[nn*(ny*(2*xoff-1-i)+j)+n];
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_X)
    /* Periodic. avoid communication to myself */
    for (n=0;n<nn;n++){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++){
	  f[n][nx*j+(nx-1-i)]=f[n][nx*j+(2*xoff-1-i)];
	  f[n][nx*j+i]=f[n][nx*j+(nx-2*xoff+i)];
	}
      }
    }
#endif
  }

  /* YBC */
#if (PRDC_Y)
  rankl=((mpi_rank / mpi_numx) == 0)?(mpi_rank+mpi_numx*(mpi_numy-1)):(mpi_rank-mpi_numx);
  rankh=((mpi_rank / mpi_numx) == (mpi_numy-1))?(mpi_rank-mpi_numx*(mpi_numy-1)):(mpi_rank+mpi_numx);
#else
  rankl=((mpi_rank / mpi_numx) == 0)?(MPI_PROC_NULL):(mpi_rank-mpi_numx);
  rankh=((mpi_rank / mpi_numx) == (mpi_numy-1))?(MPI_PROC_NULL):(mpi_rank+mpi_numx);
#endif
  if (mpi_numy != 1){
    ntot=nn*nx*yoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	for (n=0;n<nn;n++){
	  fold[nn*(nx*j+i)+n]=f[n][nx*(yoff+j)+i];
	  fold[nn*(nx*(2*yoff-1-j)+i)+n]=f[n][nx*(ny-yoff-1-j)+i];
	  fcpy[nn*(nx*j+i)+n]=f[n][nx*j+i];
	  fcpy[nn*(nx*(2*yoff-1-j)+i)+n]=f[n][nx*(ny-1-j)+i];
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	for (n=0;n<nn;n++){
	  f[n][nx*j+i]=fcpy[nn*(nx*j+i)+n];
	  f[n][nx*(ny-1-j)+i]=fcpy[nn*(nx*(2*yoff-1-j)+i)+n];
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_Y)
    /* Periodic. avoid communication to myself */
    for (n=0;n<nn;n++){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  f[n][nx*(ny-1-j)+i]=f[n][nx*(2*yoff-1-j)+i];
	  f[n][nx*j+i]=f[n][nx*(ny-2*yoff+j)+i];
	}
      }
    }
#endif
  }

}

void mpi_xbc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D X Dirichlet or Neumann BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1) or Neumann (+1). if dn==0, nothing to do */
{
  if (dn != 0 && PRDC_X == 0){
    int i,j;
    /* Left */
    if ((mpi_rank % mpi_numx) == 0){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++) f[nx*j+i]=dn*f[nx*j+(2*xoff-1+st)-i];
      }
    }
    /* Right */
    if ((mpi_rank % mpi_numx) == (mpi_numx-1)){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff-st;i++) f[nx*j+(nx-1-i)]=dn*f[nx*j+(nx-2*xoff+st)+i];
      }
    }
  }
}

void mpi_ybc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D Y Dirichlet or Neumann BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1) or Neumann (+1). if dn==0, nothing to do */
{
  if (dn != 0 && PRDC_Y == 0){
    int i,j;
    /* Left */
    if (mpi_rank/mpi_numx == 0){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++) f[nx*j+i]=dn*f[nx*((2*yoff-1+st)-j)+i];
      }
    }
    /* Right */
    if (mpi_rank/mpi_numx == (mpi_numy-1)){
      for (j=0;j<yoff-st;j++){
	for (i=0;i<nx;i++) f[nx*(ny-1-j)+i]=dn*f[nx*((ny-2*yoff+st)+j)+i];
      }
    }
  }
}
