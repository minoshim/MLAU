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

#define PRDC_X (1)
#define PRDC_Y (0)

void mpi_sdrv2d(double *f, int nx, int ny, int xoff, int yoff, 
		int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D MPI SendRecv */
{
  int i,j;
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
    ntot=ny*xoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (i=0;i<xoff;i++){
      for (j=0;j<ny;j++){	/* Transpose */
	fold[ny*i+j]=f[nx*j+(xoff+i)];
	fold[ny*(2*xoff-1-i)+j]=f[nx*j+(nx-xoff-1-i)];
	fcpy[ny*i+j]=f[nx*j+i];
	fcpy[ny*(2*xoff-1-i)+j]=f[nx*j+(nx-1-i)];
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
	f[nx*j+i]=fcpy[ny*i+j];
	f[nx*j+(nx-1-i)]=fcpy[ny*(2*xoff-1-i)+j];
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_X)
    /* Periodic. avoid communication to myself */
    for (j=0;j<ny;j++){
      for (i=0;i<xoff;i++){
	f[nx*j+(nx-1-i)]=f[nx*j+(2*xoff-1-i)];
	f[nx*j+i]=f[nx*j+(nx-2*xoff+i)];
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
    ntot=nx*yoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	fold[nx*j+i]=f[nx*(yoff+j)+i];
	fold[nx*(2*yoff-1-j)+i]=f[nx*(ny-yoff-1-j)+i];
	fcpy[nx*j+i]=f[nx*j+i];
	fcpy[nx*(2*yoff-1-j)+i]=f[nx*(ny-1-j)+i];
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
	f[nx*j+i]=fcpy[nx*j+i];
	f[nx*(ny-1-j)+i]=fcpy[nx*(2*yoff-1-j)+i];
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_Y)
    /* Periodic. avoid communication to myself */
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	f[nx*(ny-1-j)+i]=f[nx*(2*yoff-1-j)+i];
	f[nx*j+i]=f[nx*(ny-2*yoff+j)+i];
      }
    }
#endif
  }

}
