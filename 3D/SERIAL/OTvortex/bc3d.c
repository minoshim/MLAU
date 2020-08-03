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

void bc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff,
	  int stx, int dnx, int sty, int dny, int stz, int dnz)
/* 3D Boundary condition */

/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center), else 0 */
/* dn: 0 for periodic, -1 for Dirichlet, +1 for Neumann condition */
{
  int i,j,k;

  if (dnx == 0){
    /* Periodic */
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++){
	  f[nx*(ny*k+j)+(nx-1-i)]=f[nx*(ny*k+j)+( 2*xoff-1-i)];
	  f[nx*(ny*k+j)+(     i)]=f[nx*(ny*k+j)+(nx-2*xoff+i)];
	}
      }
    }
  } else{
    /* Dirichlet (dn = -1) or Neumann (dn = +1) */
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++){
	  /* Left */
	  f[nx*(ny*k+j)+(     i)]=dnx*f[nx*(ny*k+j)+( 2*xoff-1+stx-i)];
	}
	for (i=0;i<xoff-stx;i++){
	  /* Right */
	  f[nx*(ny*k+j)+(nx-1-i)]=dnx*f[nx*(ny*k+j)+(nx-2*xoff+stx+i)];
	}
      }
    }
  }
  
  if (dny == 0){
    /* Periodic */
    for (k=0;k<nz;k++){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  f[nx*(ny*k+(ny-1-j))+i]=f[nx*(ny*k+( 2*yoff-1-j))+i];
	  f[nx*(ny*k+(     j))+i]=f[nx*(ny*k+(ny-2*yoff+j))+i];
	}
      }
    }
  } else{
    /* Dirichlet (dn = -1) or Neumann (dn = +1) */
    for (k=0;k<nz;k++){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  /* Left */
	  f[nx*(ny*k+(     j))+i]=dny*f[nx*(ny*k+( 2*yoff-1+sty-j))+i];	
	}
      }
      for (j=0;j<yoff-sty;j++){
	for (i=0;i<nx;i++){
	  /* Right */
	  f[nx*(ny*k+(ny-1-j))+i]=dny*f[nx*(ny*k+(ny-2*yoff+sty+j))+i];
	}
      }
    }
  }

  if (dnz == 0){
    /* Periodic */
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  f[nx*(ny*(nz-1-k)+j)+i]=f[nx*(ny*( 2*zoff-1-k)+j)+i];
	  f[nx*(ny*(     k)+j)+i]=f[nx*(ny*(nz-2*zoff+k)+j)+i];
	}
      }
    }
  } else{
    /* Dirichlet (dn = -1) or Neumann (dn = +1) */
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  /* Left */
	  f[nx*(ny*(     k)+j)+i]=dnz*f[nx*(ny*( 2*zoff-1+stz-k)+j)+i];
	}
      }
    }
    for (k=0;k<zoff-stz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  /* Right */
	  f[nx*(ny*(nz-1-k)+j)+i]=dnz*f[nx*(ny*(nz-2*zoff+stz+k)+j)+i];	  
	}
      }
    }
  }

}
