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

void mpi_sdrv3d(double *f[], int nn, int nx, int ny, int nz, int xoff, int yoff, int zoff,
		int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D MPI SendRecv */
{
  int i,j,k,n;
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
    ntot=nn*ny*nz*xoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (i=0;i<xoff;i++){
      for (k=0;k<nz;k++){
	for (j=0;j<ny;j++){	/* Transpose */
	  for (n=0;n<nn;n++){
	    fold[nn*(ny*(nz*i+k)+j)+n]=f[n][nx*(ny*k+j)+(xoff+i)];
	    fold[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+n]=f[n][nx*(ny*k+j)+(nx-xoff-1-i)];
	    fcpy[nn*(ny*(nz*i+k)+j)+n]=f[n][nx*(ny*k+j)+i];
	    fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+n]=f[n][nx*(ny*k+j)+(nx-1-i)];
	  }
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
	  for (n=0;n<nn;n++){
	    f[n][nx*(ny*k+j)+i]=fcpy[nn*(ny*(nz*i+k)+j)+n];
	    f[n][nx*(ny*k+j)+(nx-1-i)]=fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+n];
	  }
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_X)
    /* Periodic. avoid communication to myself */
    for (n=0;n<nn;n++){
      for (k=0;k<nz;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<xoff;i++){
	    f[n][nx*(ny*k+j)+(nx-1-i)]=f[n][nx*(ny*k+j)+(2*xoff-1-i)];
	    f[n][nx*(ny*k+j)+i]=f[n][nx*(ny*k+j)+(nx-2*xoff+i)];
	  }
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
	  for (n=0;n<nn;n++){
	    fold[nn*(nz*(nx*j+i)+k)+n]=f[n][nx*(ny*k+(yoff+j))+i];
	    fold[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+n]=f[n][nx*(ny*k+(ny-yoff-1-j))+i];
	    fcpy[nn*(nz*(nx*j+i)+k)+n]=f[n][nx*(ny*k+j)+i];
	    fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+n]=f[n][nx*(ny*k+(ny-1-j))+i];
	  }
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
	  for (n=0;n<nn;n++){
	    f[n][nx*(ny*k+j)+i]=fcpy[nn*(nz*(nx*j+i)+k)+n];
	    f[n][nx*(ny*k+(ny-1-j))+i]=fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+n];
	  }
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_Y)
    /* Periodic. avoid communication to myself */
    for (n=0;n<nn;n++){
      for (k=0;k<nz;k++){
	for (j=0;j<yoff;j++){
	  for (i=0;i<nx;i++){
	  f[n][nx*(ny*k+(ny-1-j))+i]=f[n][nx*(ny*k+(2*yoff-1-j))+i];
	  f[n][nx*(ny*k+j)+i]=f[n][nx*(ny*k+(ny-2*yoff+j))+i];
	  }
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
	  for (n=0;n<nn;n++){
	    fold[nn*(nx*(ny*k+j)+i)+n]=f[n][nx*(ny*(zoff+k)+j)+i];
	    fold[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+n]=f[n][nx*(ny*(nz-zoff-1-k)+j)+i];
	    fcpy[nn*(nx*(ny*k+j)+i)+n]=f[n][nx*(ny*k+j)+i];
	    fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+n]=f[n][nx*(ny*(nz-1-k)+j)+i];
	  }
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
	  for (n=0;n<nn;n++){
	    f[n][nx*(ny*k+j)+i]=fcpy[nn*(nx*(ny*k+j)+i)+n];
	    f[n][nx*(ny*(nz-1-k)+j)+i]=fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+n];
	  }
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
#if (PRDC_Z)
    /* Periodic. avoid communication to myself */
    for (n=0;n<nn;n++){
      for (k=0;k<zoff;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<nx;i++){
	  f[n][nx*(ny*(nz-1-k)+j)+i]=f[n][nx*(ny*(2*zoff-1-k)+j)+i];
	  f[n][nx*(ny*k+j)+i]=f[n][nx*(ny*(nz-2*zoff+k)+j)+i];
	  }
	}
      }
    }
#endif
  }

}
