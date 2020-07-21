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

void bc2d(double *f, int nx, int ny, int xoff, int yoff,
	  int stx, int dnx, int sty, int dny)

/* 2D Boundary condition */

/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center), else 0 */
/* dn: 0 for periodic, -1 for Dirichlet, +1 for Neumann condition */
{
  int i,j;

  if (dnx == 0){
    /* Periodic */
    for (j=0;j<ny;j++){
      for (i=0;i<xoff;i++){
	f[nx*j+(nx-1-i)]=f[nx*j+(2*xoff-1-i)];
	f[nx*j+i]=f[nx*j+(nx-2*xoff+i)];
      }
    }
  } else{
    /* Dirichlet (dn = -1) or Neumann (dn = +1) */
    for (j=0;j<ny;j++){
      for (i=0;i<xoff;i++){
	f[nx*j+i]=dnx*f[nx*j+(2*xoff-1+stx)-i];
	f[nx*j+(nx-1-i)]=dnx*f[nx*j+(nx-2*xoff+stx)+i];
      }
    }
  }
  
  if (dny == 0){
    /* Periodic */
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	f[nx*(ny-1-j)+i]=f[nx*(2*yoff-1-j)+i];
	f[nx*j+i]=f[nx*(ny-2*yoff+j)+i];
      }
    }
  } else{
    /* Dirichlet (dn = -1) or Neumann (dn = +1) */
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	f[nx*j+i]=dny*f[nx*(2*yoff-1+sty-j)+i];	
	f[nx*(ny-1-j)+i]=dny*f[nx*(ny-2*yoff+sty+j)+i];
      }
    }
  }

}
