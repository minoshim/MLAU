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
#define PRDC_Y (1)
#define PRDC_Z (1)

void mpi_sdrv3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff,
		int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);
void mpi_sdrv3d_08(double *f0, double *f1, double *f2, double *f3, 
		   double *f4, double *f5, double *f6, double *f7,
		   int nx, int ny, int nz, int xoff, int yoff, int zoff,
		   int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);

void mpi_sdrv3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff,
		int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D MPI SendRecv */
{
  unsigned long i,j,k;
  int m_xy=mpi_numx*mpi_numy;
  int mpi_tag=0;
  int rankl,rankh;
  int ntot,ntot2;
  MPI_Status r_stat;
  double *fold,*fcpy;

  /* XBC */
#if (PRDC_X)
  rankl=(((mpi_rank%m_xy)%mpi_numx) == 0)?(mpi_rank+(mpi_numx-1)):(mpi_rank-1);
  rankh=(((mpi_rank%m_xy)%mpi_numx) == (mpi_numx-1))?(mpi_rank-(mpi_numx-1)):(mpi_rank+1);
#else
  rankl=(((mpi_rank%m_xy)%mpi_numx) == 0)?(MPI_PROC_NULL):(mpi_rank-1);
  rankh=(((mpi_rank%m_xy)%mpi_numx) == (mpi_numx-1))?(MPI_PROC_NULL):(mpi_rank+1);
#endif  
  if (mpi_numx != 1){
    ntot=ny*nz*xoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (i=0;i<xoff;i++){
      for (k=0;k<nz;k++){
	for (j=0;j<ny;j++){	/* Transpose */
	  fold[ny*(nz*i+k)+j]=f[nx*(ny*k+j)+(xoff+i)];
	  fold[ny*(nz*(2*xoff-1-i)+k)+j]=f[nx*(ny*k+j)+(nx-xoff-1-i)];
	  fcpy[ny*(nz*i+k)+j]=f[nx*(ny*k+j)+i];
	  fcpy[ny*(nz*(2*xoff-1-i)+k)+j]=f[nx*(ny*k+j)+(nx-1-i)];
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++){
	  f[nx*(ny*k+j)+i]=fcpy[ny*(nz*i+k)+j];
	  f[nx*(ny*k+j)+(nx-1-i)]=fcpy[ny*(nz*(2*xoff-1-i)+k)+j];
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_X)
    /* Periodic. avoid communication to myself */
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++){
	  f[nx*(ny*k+j)+(nx-1-i)]=f[nx*(ny*k+j)+(2*xoff-1-i)];
	  f[nx*(ny*k+j)+i]=f[nx*(ny*k+j)+(nx-2*xoff+i)];
	}
      }
    }
#endif
  }
  
  /* YBC */
#if (PRDC_Y)
  rankl=(((mpi_rank%m_xy)/mpi_numx) == 0)?(mpi_rank+mpi_numx*(mpi_numy-1)):(mpi_rank-mpi_numx);
  rankh=(((mpi_rank%m_xy)/mpi_numx) == (mpi_numy-1))?(mpi_rank-mpi_numx*(mpi_numy-1)):(mpi_rank+mpi_numx);
#else
  rankl=(((mpi_rank%m_xy)/mpi_numx) == 0)?(MPI_PROC_NULL):(mpi_rank-mpi_numx);
  rankh=(((mpi_rank%m_xy)/mpi_numx) == (mpi_numy-1))?(MPI_PROC_NULL):(mpi_rank+mpi_numx);
#endif
  if (mpi_numy != 1){
    ntot=nz*nx*yoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	for (k=0;k<nz;k++){	/* Transpose */
	  fold[nz*(nx*j+i)+k]=f[nx*(ny*k+(yoff+j))+i];
	  fold[nz*(nx*(2*yoff-1-j)+i)+k]=f[nx*(ny*k+(ny-yoff-1-j))+i];
	  fcpy[nz*(nx*j+i)+k]=f[nx*(ny*k+j)+i];
	  fcpy[nz*(nx*(2*yoff-1-j)+i)+k]=f[nx*(ny*k+(ny-1-j))+i];
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (k=0;k<nz;k++){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  f[nx*(ny*k+j)+i]=fcpy[nz*(nx*j+i)+k];
	  f[nx*(ny*k+(ny-1-j))+i]=fcpy[nz*(nx*(2*yoff-1-j)+i)+k];
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_Y)
    /* Periodic. avoid communication to myself */
    for (k=0;k<nz;k++){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  f[nx*(ny*k+(ny-1-j))+i]=f[nx*(ny*k+(2*yoff-1-j))+i];
	  f[nx*(ny*k+j)+i]=f[nx*(ny*k+(ny-2*yoff+j))+i];
	}
      }
    }
#endif
  }

  /* ZBC */
#if (PRDC_Z)
  rankl=((mpi_rank/m_xy) == 0)?(mpi_rank+m_xy*(mpi_numz-1)):(mpi_rank-m_xy);
  rankh=((mpi_rank/m_xy) == (mpi_numz-1))?(mpi_rank-m_xy*(mpi_numz-1)):(mpi_rank+m_xy);
#else
  rankl=((mpi_rank/m_xy) == 0)?(MPI_PROC_NULL):(mpi_rank-m_xy);
  rankh=((mpi_rank/m_xy) == (mpi_numz-1))?(MPI_PROC_NULL):(mpi_rank+m_xy);
#endif
  if (mpi_numz != 1){
    ntot=nx*ny*zoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  fold[nx*(ny*k+j)+i]=f[nx*(ny*(zoff+k)+j)+i];
	  fold[nx*(ny*(2*zoff-1-k)+j)+i]=f[nx*(ny*(nz-zoff-1-k)+j)+i];
	  fcpy[nx*(ny*k+j)+i]=f[nx*(ny*k+j)+i];
	  fcpy[nx*(ny*(2*zoff-1-k)+j)+i]=f[nx*(ny*(nz-1-k)+j)+i];
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  f[nx*(ny*k+j)+i]=fcpy[nx*(ny*k+j)+i];
	  f[nx*(ny*(nz-1-k)+j)+i]=fcpy[nx*(ny*(2*zoff-1-k)+j)+i];
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_Z)
    /* Periodic. avoid communication to myself */
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  f[nx*(ny*(nz-1-k)+j)+i]=f[nx*(ny*(2*zoff-1-k)+j)+i];
	  f[nx*(ny*k+j)+i]=f[nx*(ny*(nz-2*zoff+k)+j)+i];
	}
      }
    }
#endif
  }

}

void mpi_sdrv3d_08(double *f0, double *f1, double *f2, double *f3, 
		   double *f4, double *f5, double *f6, double *f7,
		   int nx, int ny, int nz, int xoff, int yoff, int zoff,
		   int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D MPI SendRecv */
/* Pack 8 variables */
{
  unsigned long i,j,k;
  int m_xy=mpi_numx*mpi_numy;
  int nn=8;
  int mpi_tag=0;
  int rankl,rankh;
  int ntot,ntot2;
  MPI_Status r_stat;
  double *fold,*fcpy;

  /* XBC */
#if (PRDC_X)
  rankl=(((mpi_rank%m_xy)%mpi_numx) == 0)?(mpi_rank+(mpi_numx-1)):(mpi_rank-1);
  rankh=(((mpi_rank%m_xy)%mpi_numx) == (mpi_numx-1))?(mpi_rank-(mpi_numx-1)):(mpi_rank+1);
#else
  rankl=(((mpi_rank%m_xy)%mpi_numx) == 0)?(MPI_PROC_NULL):(mpi_rank-1);
  rankh=(((mpi_rank%m_xy)%mpi_numx) == (mpi_numx-1))?(MPI_PROC_NULL):(mpi_rank+1);
#endif  
  if (mpi_numx != 1){
    ntot=nn*ny*nz*xoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (i=0;i<xoff;i++){
      for (k=0;k<nz;k++){
	for (j=0;j<ny;j++){	/* Transpose */
	  fold[nn*(ny*(nz*i+k)+j)+0]=f0[nx*(ny*k+j)+(xoff+i)];
	  fold[nn*(ny*(nz*i+k)+j)+1]=f1[nx*(ny*k+j)+(xoff+i)];
	  fold[nn*(ny*(nz*i+k)+j)+2]=f2[nx*(ny*k+j)+(xoff+i)];
	  fold[nn*(ny*(nz*i+k)+j)+3]=f3[nx*(ny*k+j)+(xoff+i)];
	  fold[nn*(ny*(nz*i+k)+j)+4]=f4[nx*(ny*k+j)+(xoff+i)];
	  fold[nn*(ny*(nz*i+k)+j)+5]=f5[nx*(ny*k+j)+(xoff+i)];
	  fold[nn*(ny*(nz*i+k)+j)+6]=f6[nx*(ny*k+j)+(xoff+i)];
	  fold[nn*(ny*(nz*i+k)+j)+7]=f7[nx*(ny*k+j)+(xoff+i)];
	  fold[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+0]=f0[nx*(ny*k+j)+(nx-xoff-1-i)];
	  fold[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+1]=f1[nx*(ny*k+j)+(nx-xoff-1-i)];
	  fold[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+2]=f2[nx*(ny*k+j)+(nx-xoff-1-i)];
	  fold[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+3]=f3[nx*(ny*k+j)+(nx-xoff-1-i)];
	  fold[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+4]=f4[nx*(ny*k+j)+(nx-xoff-1-i)];
	  fold[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+5]=f5[nx*(ny*k+j)+(nx-xoff-1-i)];
	  fold[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+6]=f6[nx*(ny*k+j)+(nx-xoff-1-i)];
	  fold[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+7]=f7[nx*(ny*k+j)+(nx-xoff-1-i)];
	  fcpy[nn*(ny*(nz*i+k)+j)+0]=f0[nx*(ny*k+j)+i];
	  fcpy[nn*(ny*(nz*i+k)+j)+1]=f1[nx*(ny*k+j)+i];
	  fcpy[nn*(ny*(nz*i+k)+j)+2]=f2[nx*(ny*k+j)+i];
	  fcpy[nn*(ny*(nz*i+k)+j)+3]=f3[nx*(ny*k+j)+i];
	  fcpy[nn*(ny*(nz*i+k)+j)+4]=f4[nx*(ny*k+j)+i];
	  fcpy[nn*(ny*(nz*i+k)+j)+5]=f5[nx*(ny*k+j)+i];
	  fcpy[nn*(ny*(nz*i+k)+j)+6]=f6[nx*(ny*k+j)+i];
	  fcpy[nn*(ny*(nz*i+k)+j)+7]=f7[nx*(ny*k+j)+i];
	  fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+0]=f0[nx*(ny*k+j)+(nx-1-i)];
	  fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+1]=f1[nx*(ny*k+j)+(nx-1-i)];
	  fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+2]=f2[nx*(ny*k+j)+(nx-1-i)];
	  fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+3]=f3[nx*(ny*k+j)+(nx-1-i)];
	  fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+4]=f4[nx*(ny*k+j)+(nx-1-i)];
	  fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+5]=f5[nx*(ny*k+j)+(nx-1-i)];
	  fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+6]=f6[nx*(ny*k+j)+(nx-1-i)];
	  fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+7]=f7[nx*(ny*k+j)+(nx-1-i)];
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++){
	  f0[nx*(ny*k+j)+i]=fcpy[nn*(ny*(nz*i+k)+j)+0];
	  f1[nx*(ny*k+j)+i]=fcpy[nn*(ny*(nz*i+k)+j)+1];
	  f2[nx*(ny*k+j)+i]=fcpy[nn*(ny*(nz*i+k)+j)+2];
	  f3[nx*(ny*k+j)+i]=fcpy[nn*(ny*(nz*i+k)+j)+3];
	  f4[nx*(ny*k+j)+i]=fcpy[nn*(ny*(nz*i+k)+j)+4];
	  f5[nx*(ny*k+j)+i]=fcpy[nn*(ny*(nz*i+k)+j)+5];
	  f6[nx*(ny*k+j)+i]=fcpy[nn*(ny*(nz*i+k)+j)+6];
	  f7[nx*(ny*k+j)+i]=fcpy[nn*(ny*(nz*i+k)+j)+7];
	  f0[nx*(ny*k+j)+(nx-1-i)]=fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+0];
	  f1[nx*(ny*k+j)+(nx-1-i)]=fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+1];
	  f2[nx*(ny*k+j)+(nx-1-i)]=fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+2];
	  f3[nx*(ny*k+j)+(nx-1-i)]=fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+3];
	  f4[nx*(ny*k+j)+(nx-1-i)]=fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+4];
	  f5[nx*(ny*k+j)+(nx-1-i)]=fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+5];
	  f6[nx*(ny*k+j)+(nx-1-i)]=fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+6];
	  f7[nx*(ny*k+j)+(nx-1-i)]=fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+7];
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_X)
    /* Periodic. avoid communication to myself */
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++){
	  f0[nx*(ny*k+j)+(nx-1-i)]=f0[nx*(ny*k+j)+(2*xoff-1-i)];
	  f0[nx*(ny*k+j)+i]=f0[nx*(ny*k+j)+(nx-2*xoff+i)];
	  f1[nx*(ny*k+j)+(nx-1-i)]=f1[nx*(ny*k+j)+(2*xoff-1-i)];
	  f1[nx*(ny*k+j)+i]=f1[nx*(ny*k+j)+(nx-2*xoff+i)];
	  f2[nx*(ny*k+j)+(nx-1-i)]=f2[nx*(ny*k+j)+(2*xoff-1-i)];
	  f2[nx*(ny*k+j)+i]=f2[nx*(ny*k+j)+(nx-2*xoff+i)];
	  f3[nx*(ny*k+j)+(nx-1-i)]=f3[nx*(ny*k+j)+(2*xoff-1-i)];
	  f3[nx*(ny*k+j)+i]=f3[nx*(ny*k+j)+(nx-2*xoff+i)];
	  f4[nx*(ny*k+j)+(nx-1-i)]=f4[nx*(ny*k+j)+(2*xoff-1-i)];
	  f4[nx*(ny*k+j)+i]=f4[nx*(ny*k+j)+(nx-2*xoff+i)];
	  f5[nx*(ny*k+j)+(nx-1-i)]=f5[nx*(ny*k+j)+(2*xoff-1-i)];
	  f5[nx*(ny*k+j)+i]=f5[nx*(ny*k+j)+(nx-2*xoff+i)];
	  f6[nx*(ny*k+j)+(nx-1-i)]=f6[nx*(ny*k+j)+(2*xoff-1-i)];
	  f6[nx*(ny*k+j)+i]=f6[nx*(ny*k+j)+(nx-2*xoff+i)];
	  f7[nx*(ny*k+j)+(nx-1-i)]=f7[nx*(ny*k+j)+(2*xoff-1-i)];
	  f7[nx*(ny*k+j)+i]=f7[nx*(ny*k+j)+(nx-2*xoff+i)];
	}
      }
    }
#endif
  }
  
  /* YBC */
#if (PRDC_Y)
  rankl=(((mpi_rank%m_xy)/mpi_numx) == 0)?(mpi_rank+mpi_numx*(mpi_numy-1)):(mpi_rank-mpi_numx);
  rankh=(((mpi_rank%m_xy)/mpi_numx) == (mpi_numy-1))?(mpi_rank-mpi_numx*(mpi_numy-1)):(mpi_rank+mpi_numx);
#else
  rankl=(((mpi_rank%m_xy)/mpi_numx) == 0)?(MPI_PROC_NULL):(mpi_rank-mpi_numx);
  rankh=(((mpi_rank%m_xy)/mpi_numx) == (mpi_numy-1))?(MPI_PROC_NULL):(mpi_rank+mpi_numx);
#endif
  if (mpi_numy != 1){
    ntot=nn*nz*nx*yoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	for (k=0;k<nz;k++){	/* Transpose */
	  fold[nn*(nz*(nx*j+i)+k)+0]=f0[nx*(ny*k+(yoff+j))+i];
	  fold[nn*(nz*(nx*j+i)+k)+1]=f1[nx*(ny*k+(yoff+j))+i];
	  fold[nn*(nz*(nx*j+i)+k)+2]=f2[nx*(ny*k+(yoff+j))+i];
	  fold[nn*(nz*(nx*j+i)+k)+3]=f3[nx*(ny*k+(yoff+j))+i];
	  fold[nn*(nz*(nx*j+i)+k)+4]=f4[nx*(ny*k+(yoff+j))+i];
	  fold[nn*(nz*(nx*j+i)+k)+5]=f5[nx*(ny*k+(yoff+j))+i];
	  fold[nn*(nz*(nx*j+i)+k)+6]=f6[nx*(ny*k+(yoff+j))+i];
	  fold[nn*(nz*(nx*j+i)+k)+7]=f7[nx*(ny*k+(yoff+j))+i];
	  fold[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+0]=f0[nx*(ny*k+(ny-yoff-1-j))+i];
	  fold[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+1]=f1[nx*(ny*k+(ny-yoff-1-j))+i];
	  fold[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+2]=f2[nx*(ny*k+(ny-yoff-1-j))+i];
	  fold[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+3]=f3[nx*(ny*k+(ny-yoff-1-j))+i];
	  fold[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+4]=f4[nx*(ny*k+(ny-yoff-1-j))+i];
	  fold[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+5]=f5[nx*(ny*k+(ny-yoff-1-j))+i];
	  fold[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+6]=f6[nx*(ny*k+(ny-yoff-1-j))+i];
	  fold[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+7]=f7[nx*(ny*k+(ny-yoff-1-j))+i];
	  fcpy[nn*(nz*(nx*j+i)+k)+0]=f0[nx*(ny*k+j)+i];
	  fcpy[nn*(nz*(nx*j+i)+k)+1]=f1[nx*(ny*k+j)+i];
	  fcpy[nn*(nz*(nx*j+i)+k)+2]=f2[nx*(ny*k+j)+i];
	  fcpy[nn*(nz*(nx*j+i)+k)+3]=f3[nx*(ny*k+j)+i];
	  fcpy[nn*(nz*(nx*j+i)+k)+4]=f4[nx*(ny*k+j)+i];
	  fcpy[nn*(nz*(nx*j+i)+k)+5]=f5[nx*(ny*k+j)+i];
	  fcpy[nn*(nz*(nx*j+i)+k)+6]=f6[nx*(ny*k+j)+i];
	  fcpy[nn*(nz*(nx*j+i)+k)+7]=f7[nx*(ny*k+j)+i];
	  fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+0]=f0[nx*(ny*k+(ny-1-j))+i];
	  fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+1]=f1[nx*(ny*k+(ny-1-j))+i];
	  fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+2]=f2[nx*(ny*k+(ny-1-j))+i];
	  fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+3]=f3[nx*(ny*k+(ny-1-j))+i];
	  fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+4]=f4[nx*(ny*k+(ny-1-j))+i];
	  fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+5]=f5[nx*(ny*k+(ny-1-j))+i];
	  fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+6]=f6[nx*(ny*k+(ny-1-j))+i];
	  fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+7]=f7[nx*(ny*k+(ny-1-j))+i];
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (k=0;k<nz;k++){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  f0[nx*(ny*k+j)+i]=fcpy[nn*(nz*(nx*j+i)+k)+0];
	  f1[nx*(ny*k+j)+i]=fcpy[nn*(nz*(nx*j+i)+k)+1];
	  f2[nx*(ny*k+j)+i]=fcpy[nn*(nz*(nx*j+i)+k)+2];
	  f3[nx*(ny*k+j)+i]=fcpy[nn*(nz*(nx*j+i)+k)+3];
	  f4[nx*(ny*k+j)+i]=fcpy[nn*(nz*(nx*j+i)+k)+4];
	  f5[nx*(ny*k+j)+i]=fcpy[nn*(nz*(nx*j+i)+k)+5];
	  f6[nx*(ny*k+j)+i]=fcpy[nn*(nz*(nx*j+i)+k)+6];
	  f7[nx*(ny*k+j)+i]=fcpy[nn*(nz*(nx*j+i)+k)+7];
	  f0[nx*(ny*k+(ny-1-j))+i]=fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+0];
	  f1[nx*(ny*k+(ny-1-j))+i]=fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+1];
	  f2[nx*(ny*k+(ny-1-j))+i]=fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+2];
	  f3[nx*(ny*k+(ny-1-j))+i]=fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+3];
	  f4[nx*(ny*k+(ny-1-j))+i]=fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+4];
	  f5[nx*(ny*k+(ny-1-j))+i]=fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+5];
	  f6[nx*(ny*k+(ny-1-j))+i]=fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+6];
	  f7[nx*(ny*k+(ny-1-j))+i]=fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+7];
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_Y)
    /* Periodic. avoid communication to myself */
    for (k=0;k<nz;k++){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  f0[nx*(ny*k+(ny-1-j))+i]=f0[nx*(ny*k+(2*yoff-1-j))+i];
	  f0[nx*(ny*k+j)+i]=f0[nx*(ny*k+(ny-2*yoff+j))+i];
	  f1[nx*(ny*k+(ny-1-j))+i]=f1[nx*(ny*k+(2*yoff-1-j))+i];
	  f1[nx*(ny*k+j)+i]=f1[nx*(ny*k+(ny-2*yoff+j))+i];
	  f2[nx*(ny*k+(ny-1-j))+i]=f2[nx*(ny*k+(2*yoff-1-j))+i];
	  f2[nx*(ny*k+j)+i]=f2[nx*(ny*k+(ny-2*yoff+j))+i];
	  f3[nx*(ny*k+(ny-1-j))+i]=f3[nx*(ny*k+(2*yoff-1-j))+i];
	  f3[nx*(ny*k+j)+i]=f3[nx*(ny*k+(ny-2*yoff+j))+i];
	  f4[nx*(ny*k+(ny-1-j))+i]=f4[nx*(ny*k+(2*yoff-1-j))+i];
	  f4[nx*(ny*k+j)+i]=f4[nx*(ny*k+(ny-2*yoff+j))+i];
	  f5[nx*(ny*k+(ny-1-j))+i]=f5[nx*(ny*k+(2*yoff-1-j))+i];
	  f5[nx*(ny*k+j)+i]=f5[nx*(ny*k+(ny-2*yoff+j))+i];
	  f6[nx*(ny*k+(ny-1-j))+i]=f6[nx*(ny*k+(2*yoff-1-j))+i];
	  f6[nx*(ny*k+j)+i]=f6[nx*(ny*k+(ny-2*yoff+j))+i];
	  f7[nx*(ny*k+(ny-1-j))+i]=f7[nx*(ny*k+(2*yoff-1-j))+i];
	  f7[nx*(ny*k+j)+i]=f7[nx*(ny*k+(ny-2*yoff+j))+i];
	}
      }
    }
#endif
  }

  /* ZBC */
#if (PRDC_Z)
  rankl=((mpi_rank/m_xy) == 0)?(mpi_rank+m_xy*(mpi_numz-1)):(mpi_rank-m_xy);
  rankh=((mpi_rank/m_xy) == (mpi_numz-1))?(mpi_rank-m_xy*(mpi_numz-1)):(mpi_rank+m_xy);
#else
  rankl=((mpi_rank/m_xy) == 0)?(MPI_PROC_NULL):(mpi_rank-m_xy);
  rankh=((mpi_rank/m_xy) == (mpi_numz-1))?(MPI_PROC_NULL):(mpi_rank+m_xy);
#endif
  if (mpi_numz != 1){
    ntot=nn*nx*ny*zoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  fold[nn*(nx*(ny*k+j)+i)+0]=f0[nx*(ny*(zoff+k)+j)+i];
	  fold[nn*(nx*(ny*k+j)+i)+1]=f1[nx*(ny*(zoff+k)+j)+i];
	  fold[nn*(nx*(ny*k+j)+i)+2]=f2[nx*(ny*(zoff+k)+j)+i];
	  fold[nn*(nx*(ny*k+j)+i)+3]=f3[nx*(ny*(zoff+k)+j)+i];
	  fold[nn*(nx*(ny*k+j)+i)+4]=f4[nx*(ny*(zoff+k)+j)+i];
	  fold[nn*(nx*(ny*k+j)+i)+5]=f5[nx*(ny*(zoff+k)+j)+i];
	  fold[nn*(nx*(ny*k+j)+i)+6]=f6[nx*(ny*(zoff+k)+j)+i];
	  fold[nn*(nx*(ny*k+j)+i)+7]=f7[nx*(ny*(zoff+k)+j)+i];
	  fold[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+0]=f0[nx*(ny*(nz-zoff-1-k)+j)+i];
	  fold[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+1]=f1[nx*(ny*(nz-zoff-1-k)+j)+i];
	  fold[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+2]=f2[nx*(ny*(nz-zoff-1-k)+j)+i];
	  fold[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+3]=f3[nx*(ny*(nz-zoff-1-k)+j)+i];
	  fold[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+4]=f4[nx*(ny*(nz-zoff-1-k)+j)+i];
	  fold[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+5]=f5[nx*(ny*(nz-zoff-1-k)+j)+i];
	  fold[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+6]=f6[nx*(ny*(nz-zoff-1-k)+j)+i];
	  fold[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+7]=f7[nx*(ny*(nz-zoff-1-k)+j)+i];
	  fcpy[nn*(nx*(ny*k+j)+i)+0]=f0[nx*(ny*k+j)+i];
	  fcpy[nn*(nx*(ny*k+j)+i)+1]=f1[nx*(ny*k+j)+i];
	  fcpy[nn*(nx*(ny*k+j)+i)+2]=f2[nx*(ny*k+j)+i];
	  fcpy[nn*(nx*(ny*k+j)+i)+3]=f3[nx*(ny*k+j)+i];
	  fcpy[nn*(nx*(ny*k+j)+i)+4]=f4[nx*(ny*k+j)+i];
	  fcpy[nn*(nx*(ny*k+j)+i)+5]=f5[nx*(ny*k+j)+i];
	  fcpy[nn*(nx*(ny*k+j)+i)+6]=f6[nx*(ny*k+j)+i];
	  fcpy[nn*(nx*(ny*k+j)+i)+7]=f7[nx*(ny*k+j)+i];
	  fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+0]=f0[nx*(ny*(nz-1-k)+j)+i];
	  fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+1]=f1[nx*(ny*(nz-1-k)+j)+i];
	  fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+2]=f2[nx*(ny*(nz-1-k)+j)+i];
	  fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+3]=f3[nx*(ny*(nz-1-k)+j)+i];
	  fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+4]=f4[nx*(ny*(nz-1-k)+j)+i];
	  fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+5]=f5[nx*(ny*(nz-1-k)+j)+i];
	  fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+6]=f6[nx*(ny*(nz-1-k)+j)+i];
	  fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+7]=f7[nx*(ny*(nz-1-k)+j)+i];
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  f0[nx*(ny*k+j)+i]=fcpy[nn*(nx*(ny*k+j)+i)+0];
	  f1[nx*(ny*k+j)+i]=fcpy[nn*(nx*(ny*k+j)+i)+1];
	  f2[nx*(ny*k+j)+i]=fcpy[nn*(nx*(ny*k+j)+i)+2];
	  f3[nx*(ny*k+j)+i]=fcpy[nn*(nx*(ny*k+j)+i)+3];
	  f4[nx*(ny*k+j)+i]=fcpy[nn*(nx*(ny*k+j)+i)+4];
	  f5[nx*(ny*k+j)+i]=fcpy[nn*(nx*(ny*k+j)+i)+5];
	  f6[nx*(ny*k+j)+i]=fcpy[nn*(nx*(ny*k+j)+i)+6];
	  f7[nx*(ny*k+j)+i]=fcpy[nn*(nx*(ny*k+j)+i)+7];
	  f0[nx*(ny*(nz-1-k)+j)+i]=fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+0];
	  f1[nx*(ny*(nz-1-k)+j)+i]=fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+1];
	  f2[nx*(ny*(nz-1-k)+j)+i]=fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+2];
	  f3[nx*(ny*(nz-1-k)+j)+i]=fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+3];
	  f4[nx*(ny*(nz-1-k)+j)+i]=fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+4];
	  f5[nx*(ny*(nz-1-k)+j)+i]=fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+5];
	  f6[nx*(ny*(nz-1-k)+j)+i]=fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+6];
	  f7[nx*(ny*(nz-1-k)+j)+i]=fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+7];
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_Z)
    /* Periodic. avoid communication to myself */
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  f0[nx*(ny*(nz-1-k)+j)+i]=f0[nx*(ny*(2*zoff-1-k)+j)+i];
	  f0[nx*(ny*k+j)+i]=f0[nx*(ny*(nz-2*zoff+k)+j)+i];
	  f1[nx*(ny*(nz-1-k)+j)+i]=f1[nx*(ny*(2*zoff-1-k)+j)+i];
	  f1[nx*(ny*k+j)+i]=f1[nx*(ny*(nz-2*zoff+k)+j)+i];
	  f2[nx*(ny*(nz-1-k)+j)+i]=f2[nx*(ny*(2*zoff-1-k)+j)+i];
	  f2[nx*(ny*k+j)+i]=f2[nx*(ny*(nz-2*zoff+k)+j)+i];
	  f3[nx*(ny*(nz-1-k)+j)+i]=f3[nx*(ny*(2*zoff-1-k)+j)+i];
	  f3[nx*(ny*k+j)+i]=f3[nx*(ny*(nz-2*zoff+k)+j)+i];
	  f4[nx*(ny*(nz-1-k)+j)+i]=f4[nx*(ny*(2*zoff-1-k)+j)+i];
	  f4[nx*(ny*k+j)+i]=f4[nx*(ny*(nz-2*zoff+k)+j)+i];
	  f5[nx*(ny*(nz-1-k)+j)+i]=f5[nx*(ny*(2*zoff-1-k)+j)+i];
	  f5[nx*(ny*k+j)+i]=f5[nx*(ny*(nz-2*zoff+k)+j)+i];
	  f6[nx*(ny*(nz-1-k)+j)+i]=f6[nx*(ny*(2*zoff-1-k)+j)+i];
	  f6[nx*(ny*k+j)+i]=f6[nx*(ny*(nz-2*zoff+k)+j)+i];
	  f7[nx*(ny*(nz-1-k)+j)+i]=f7[nx*(ny*(2*zoff-1-k)+j)+i];
	  f7[nx*(ny*k+j)+i]=f7[nx*(ny*(nz-2*zoff+k)+j)+i];
	}
      }
    }
#endif
  }

}
