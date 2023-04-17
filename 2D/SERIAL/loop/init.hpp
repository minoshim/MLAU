void init_grid();
void init_plasma();

void init_grid()
// Define X and Y coordinates
{
  char filname[100];
  FILE *outfil;

  sprintf(filname,"%s/x.dat",fildir);
  outfil=fopen(filname,"w");
  for (int i=0;i<nx;i++){
    x[i]=(i-xoff+0.5)*dx-0.5*lx;
    fprintf(outfil,"%.12f\n",x[i]);
  }
  fclose(outfil);

  sprintf(filname,"%s/y.dat",fildir);
  outfil=fopen(filname,"w");
  for (int j=0;j<ny;j++){
    y[j]=(j-yoff+0.5)*dy-0.5*ly;
    fprintf(outfil,"%.12f\n",y[j]);
  }
  fclose(outfil);

}

void init_plasma()
// Set initial condition
{
  // Magnetic loop advection in high beta plasma
  int i,j;
  double *azp,*azc;
  azp=new double[nx*ny];
  azc=new double[nx*ny];
  // Vector potential
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      int ss=nx*j+i;
      double r,x0,y0;
      x0=x[i]-0.5*dx;
      y0=y[j]-0.5*dy;
      r=sqrt(x0*x0+y0*y0);
      azp[ss]=(r <= rad)?a0*(rad-r):0;
      x0=x[i];
      y0=y[j];
      r=sqrt(x0*x0+y0*y0);
      azc[ss]=(r <= rad)?a0*(rad-r):0;
    }
  }
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      int ss=nx*j+i;
      double vx,vy,vz,pr,cx,cy,cz;
      ro[ss]=ro0;
      vx=vx0;
      vy=vy0;
      vz=vz0;
      pr=pr0;

      bx[ss]=by[ss]=bz[ss]=0.0;
      cx=cy=cz=0.0;
      if (j >= 1 && j <= ny-2){
	bx[ss]=+(azp[nx*(j+1)+i]-azp[nx*j+i])*idy;
	cx=+(azc[nx*(j+1)+i]-azc[nx*(j-1)+i])*idy*0.5;
      }
      if (i >= 1 && i <= nx-2){
	by[ss]=-(azp[nx*j+(i+1)]-azp[nx*j+i])*idx;
	cy=-(azc[nx*j+(i+1)]-azc[nx*j+(i-1)])*idx*0.5;
      }
      // Use the following 4th-order discretization for bx and by when ODR=4 in mhd_fd_ct_2d.c
      // if (j >= 2 && j <= ny-3){
      // 	bx[ss]=+(27.0*(azp[nx*(j+1)+i]-azp[nx*j+i])
      // 		 -(azp[nx*(j+2)+i]-azp[nx*(j-1)+i]))/24.0*idy;
      // 	cx=+(azc[nx*(j+1)+i]-azc[nx*(j-1)+i])*idy*0.5;
      // }
      // if (i >= 2 && i <= nx-3){
      // 	by[ss]=-(27.0*(azp[nx*j+(i+1)]-azp[nx*j+i])
      // 		 -(azp[nx*j+(i+2)]-azp[nx*j+(i-1)]))/24.0*idx;
      // 	cy=-(azc[nx*j+(i+1)]-azc[nx*j+(i-1)])*idx*0.5;
      // }
      
      mx[ss]=ro[ss]*vx;
      my[ss]=ro[ss]*vy;
      mz[ss]=ro[ss]*vz;
      en[ss]=pr/(gam-1)+0.5*(ro[ss]*(vx*vx+vy*vy+vz*vz)+(cx*cx+cy*cy+cz*cz));
    }
  }

  delete[] azp;
  delete[] azc;
}
