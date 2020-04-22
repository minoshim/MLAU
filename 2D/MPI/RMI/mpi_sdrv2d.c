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
		int mpi_rank, int mpi_numx, int mpi_numy);
void mpi_sdrv2d_08(double *f0, double *f1, double *f2, double *f3, 
		   double *f4, double *f5, double *f6, double *f7,
		   int nx, int ny, int xoff, int yoff, 
		   int mpi_rank, int mpi_numx, int mpi_numy);

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

void mpi_sdrv2d_08(double *f0, double *f1, double *f2, double *f3, 
		   double *f4, double *f5, double *f6, double *f7,
		   int nx, int ny, int xoff, int yoff, 
		   int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D MPI SendRecv */
/* Pack 8 variables */
{
  int i,j;
  int nn=8;
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
	fold[nn*(ny*i+j)+0]=f0[nx*j+(xoff+i)];
	fold[nn*(ny*i+j)+1]=f1[nx*j+(xoff+i)];
	fold[nn*(ny*i+j)+2]=f2[nx*j+(xoff+i)];
	fold[nn*(ny*i+j)+3]=f3[nx*j+(xoff+i)];
	fold[nn*(ny*i+j)+4]=f4[nx*j+(xoff+i)];
	fold[nn*(ny*i+j)+5]=f5[nx*j+(xoff+i)];
	fold[nn*(ny*i+j)+6]=f6[nx*j+(xoff+i)];
	fold[nn*(ny*i+j)+7]=f7[nx*j+(xoff+i)];
	fold[nn*(ny*(2*xoff-1-i)+j)+0]=f0[nx*j+(nx-xoff-1-i)];
	fold[nn*(ny*(2*xoff-1-i)+j)+1]=f1[nx*j+(nx-xoff-1-i)];
	fold[nn*(ny*(2*xoff-1-i)+j)+2]=f2[nx*j+(nx-xoff-1-i)];
	fold[nn*(ny*(2*xoff-1-i)+j)+3]=f3[nx*j+(nx-xoff-1-i)];
	fold[nn*(ny*(2*xoff-1-i)+j)+4]=f4[nx*j+(nx-xoff-1-i)];
	fold[nn*(ny*(2*xoff-1-i)+j)+5]=f5[nx*j+(nx-xoff-1-i)];
	fold[nn*(ny*(2*xoff-1-i)+j)+6]=f6[nx*j+(nx-xoff-1-i)];
	fold[nn*(ny*(2*xoff-1-i)+j)+7]=f7[nx*j+(nx-xoff-1-i)];
	fcpy[nn*(ny*i+j)+0]=f0[nx*j+i];
	fcpy[nn*(ny*i+j)+1]=f1[nx*j+i];
	fcpy[nn*(ny*i+j)+2]=f2[nx*j+i];
	fcpy[nn*(ny*i+j)+3]=f3[nx*j+i];
	fcpy[nn*(ny*i+j)+4]=f4[nx*j+i];
	fcpy[nn*(ny*i+j)+5]=f5[nx*j+i];
	fcpy[nn*(ny*i+j)+6]=f6[nx*j+i];
	fcpy[nn*(ny*i+j)+7]=f7[nx*j+i];
	fcpy[nn*(ny*(2*xoff-1-i)+j)+0]=f0[nx*j+(nx-1-i)];
	fcpy[nn*(ny*(2*xoff-1-i)+j)+1]=f1[nx*j+(nx-1-i)];
	fcpy[nn*(ny*(2*xoff-1-i)+j)+2]=f2[nx*j+(nx-1-i)];
	fcpy[nn*(ny*(2*xoff-1-i)+j)+3]=f3[nx*j+(nx-1-i)];
	fcpy[nn*(ny*(2*xoff-1-i)+j)+4]=f4[nx*j+(nx-1-i)];
	fcpy[nn*(ny*(2*xoff-1-i)+j)+5]=f5[nx*j+(nx-1-i)];
	fcpy[nn*(ny*(2*xoff-1-i)+j)+6]=f6[nx*j+(nx-1-i)];
	fcpy[nn*(ny*(2*xoff-1-i)+j)+7]=f7[nx*j+(nx-1-i)];
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
	f0[nx*j+i]=fcpy[nn*(ny*i+j)+0];
	f1[nx*j+i]=fcpy[nn*(ny*i+j)+1];
	f2[nx*j+i]=fcpy[nn*(ny*i+j)+2];
	f3[nx*j+i]=fcpy[nn*(ny*i+j)+3];
	f4[nx*j+i]=fcpy[nn*(ny*i+j)+4];
	f5[nx*j+i]=fcpy[nn*(ny*i+j)+5];
	f6[nx*j+i]=fcpy[nn*(ny*i+j)+6];
	f7[nx*j+i]=fcpy[nn*(ny*i+j)+7];
	f0[nx*j+(nx-1-i)]=fcpy[nn*(ny*(2*xoff-1-i)+j)+0];
	f1[nx*j+(nx-1-i)]=fcpy[nn*(ny*(2*xoff-1-i)+j)+1];
	f2[nx*j+(nx-1-i)]=fcpy[nn*(ny*(2*xoff-1-i)+j)+2];
	f3[nx*j+(nx-1-i)]=fcpy[nn*(ny*(2*xoff-1-i)+j)+3];
	f4[nx*j+(nx-1-i)]=fcpy[nn*(ny*(2*xoff-1-i)+j)+4];
	f5[nx*j+(nx-1-i)]=fcpy[nn*(ny*(2*xoff-1-i)+j)+5];
	f6[nx*j+(nx-1-i)]=fcpy[nn*(ny*(2*xoff-1-i)+j)+6];
	f7[nx*j+(nx-1-i)]=fcpy[nn*(ny*(2*xoff-1-i)+j)+7];
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_X)
    /* Periodic. avoid communication to myself */
    for (j=0;j<ny;j++){
      for (i=0;i<xoff;i++){
	f0[nx*j+(nx-1-i)]=f0[nx*j+(2*xoff-1-i)];
	f0[nx*j+i]=f0[nx*j+(nx-2*xoff+i)];
	f1[nx*j+(nx-1-i)]=f1[nx*j+(2*xoff-1-i)];
	f1[nx*j+i]=f1[nx*j+(nx-2*xoff+i)];
	f2[nx*j+(nx-1-i)]=f2[nx*j+(2*xoff-1-i)];
	f2[nx*j+i]=f2[nx*j+(nx-2*xoff+i)];
	f3[nx*j+(nx-1-i)]=f3[nx*j+(2*xoff-1-i)];
	f3[nx*j+i]=f3[nx*j+(nx-2*xoff+i)];
	f4[nx*j+(nx-1-i)]=f4[nx*j+(2*xoff-1-i)];
	f4[nx*j+i]=f4[nx*j+(nx-2*xoff+i)];
	f5[nx*j+(nx-1-i)]=f5[nx*j+(2*xoff-1-i)];
	f5[nx*j+i]=f5[nx*j+(nx-2*xoff+i)];
	f6[nx*j+(nx-1-i)]=f6[nx*j+(2*xoff-1-i)];
	f6[nx*j+i]=f6[nx*j+(nx-2*xoff+i)];
	f7[nx*j+(nx-1-i)]=f7[nx*j+(2*xoff-1-i)];
	f7[nx*j+i]=f7[nx*j+(nx-2*xoff+i)];
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
	fold[nn*(nx*j+i)+0]=f0[nx*(yoff+j)+i];
	fold[nn*(nx*j+i)+1]=f1[nx*(yoff+j)+i];
	fold[nn*(nx*j+i)+2]=f2[nx*(yoff+j)+i];
	fold[nn*(nx*j+i)+3]=f3[nx*(yoff+j)+i];
	fold[nn*(nx*j+i)+4]=f4[nx*(yoff+j)+i];
	fold[nn*(nx*j+i)+5]=f5[nx*(yoff+j)+i];
	fold[nn*(nx*j+i)+6]=f6[nx*(yoff+j)+i];
	fold[nn*(nx*j+i)+7]=f7[nx*(yoff+j)+i];
	fold[nn*(nx*(2*yoff-1-j)+i)+0]=f0[nx*(ny-yoff-1-j)+i];
	fold[nn*(nx*(2*yoff-1-j)+i)+1]=f1[nx*(ny-yoff-1-j)+i];
	fold[nn*(nx*(2*yoff-1-j)+i)+2]=f2[nx*(ny-yoff-1-j)+i];
	fold[nn*(nx*(2*yoff-1-j)+i)+3]=f3[nx*(ny-yoff-1-j)+i];
	fold[nn*(nx*(2*yoff-1-j)+i)+4]=f4[nx*(ny-yoff-1-j)+i];
	fold[nn*(nx*(2*yoff-1-j)+i)+5]=f5[nx*(ny-yoff-1-j)+i];
	fold[nn*(nx*(2*yoff-1-j)+i)+6]=f6[nx*(ny-yoff-1-j)+i];
	fold[nn*(nx*(2*yoff-1-j)+i)+7]=f7[nx*(ny-yoff-1-j)+i];
	fcpy[nn*(nx*j+i)+0]=f0[nx*j+i];
	fcpy[nn*(nx*j+i)+1]=f1[nx*j+i];
	fcpy[nn*(nx*j+i)+2]=f2[nx*j+i];
	fcpy[nn*(nx*j+i)+3]=f3[nx*j+i];
	fcpy[nn*(nx*j+i)+4]=f4[nx*j+i];
	fcpy[nn*(nx*j+i)+5]=f5[nx*j+i];
	fcpy[nn*(nx*j+i)+6]=f6[nx*j+i];
	fcpy[nn*(nx*j+i)+7]=f7[nx*j+i];
	fcpy[nn*(nx*(2*yoff-1-j)+i)+0]=f0[nx*(ny-1-j)+i];
	fcpy[nn*(nx*(2*yoff-1-j)+i)+1]=f1[nx*(ny-1-j)+i];
	fcpy[nn*(nx*(2*yoff-1-j)+i)+2]=f2[nx*(ny-1-j)+i];
	fcpy[nn*(nx*(2*yoff-1-j)+i)+3]=f3[nx*(ny-1-j)+i];
	fcpy[nn*(nx*(2*yoff-1-j)+i)+4]=f4[nx*(ny-1-j)+i];
	fcpy[nn*(nx*(2*yoff-1-j)+i)+5]=f5[nx*(ny-1-j)+i];
	fcpy[nn*(nx*(2*yoff-1-j)+i)+6]=f6[nx*(ny-1-j)+i];
	fcpy[nn*(nx*(2*yoff-1-j)+i)+7]=f7[nx*(ny-1-j)+i];
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
	f0[nx*j+i]=fcpy[nn*(nx*j+i)+0];
	f1[nx*j+i]=fcpy[nn*(nx*j+i)+1];
	f2[nx*j+i]=fcpy[nn*(nx*j+i)+2];
	f3[nx*j+i]=fcpy[nn*(nx*j+i)+3];
	f4[nx*j+i]=fcpy[nn*(nx*j+i)+4];
	f5[nx*j+i]=fcpy[nn*(nx*j+i)+5];
	f6[nx*j+i]=fcpy[nn*(nx*j+i)+6];
	f7[nx*j+i]=fcpy[nn*(nx*j+i)+7];
	f0[nx*(ny-1-j)+i]=fcpy[nn*(nx*(2*yoff-1-j)+i)+0];
	f1[nx*(ny-1-j)+i]=fcpy[nn*(nx*(2*yoff-1-j)+i)+1];
	f2[nx*(ny-1-j)+i]=fcpy[nn*(nx*(2*yoff-1-j)+i)+2];
	f3[nx*(ny-1-j)+i]=fcpy[nn*(nx*(2*yoff-1-j)+i)+3];
	f4[nx*(ny-1-j)+i]=fcpy[nn*(nx*(2*yoff-1-j)+i)+4];
	f5[nx*(ny-1-j)+i]=fcpy[nn*(nx*(2*yoff-1-j)+i)+5];
	f6[nx*(ny-1-j)+i]=fcpy[nn*(nx*(2*yoff-1-j)+i)+6];
	f7[nx*(ny-1-j)+i]=fcpy[nn*(nx*(2*yoff-1-j)+i)+7];
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_Y)
    /* Periodic. avoid communication to myself */
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	f0[nx*(ny-1-j)+i]=f0[nx*(2*yoff-1-j)+i];
	f0[nx*j+i]=f0[nx*(ny-2*yoff+j)+i];
	f1[nx*(ny-1-j)+i]=f1[nx*(2*yoff-1-j)+i];
	f1[nx*j+i]=f1[nx*(ny-2*yoff+j)+i];
	f2[nx*(ny-1-j)+i]=f2[nx*(2*yoff-1-j)+i];
	f2[nx*j+i]=f2[nx*(ny-2*yoff+j)+i];
	f3[nx*(ny-1-j)+i]=f3[nx*(2*yoff-1-j)+i];
	f3[nx*j+i]=f3[nx*(ny-2*yoff+j)+i];
	f4[nx*(ny-1-j)+i]=f4[nx*(2*yoff-1-j)+i];
	f4[nx*j+i]=f4[nx*(ny-2*yoff+j)+i];
	f5[nx*(ny-1-j)+i]=f5[nx*(2*yoff-1-j)+i];
	f5[nx*j+i]=f5[nx*(ny-2*yoff+j)+i];
	f6[nx*(ny-1-j)+i]=f6[nx*(2*yoff-1-j)+i];
	f6[nx*j+i]=f6[nx*(ny-2*yoff+j)+i];
	f7[nx*(ny-1-j)+i]=f7[nx*(2*yoff-1-j)+i];
	f7[nx*j+i]=f7[nx*(ny-2*yoff+j)+i];
      }
    }
#endif
  }

}
