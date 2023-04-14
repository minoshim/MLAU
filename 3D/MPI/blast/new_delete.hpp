void new_delete()
{
  static int nflg=0;
  if (nflg == 0){
    nflg=1;
    x=new double[nx];
    y=new double[ny];
    z=new double[nz];
    ro=new double[nx*ny*nz];
    mx=new double[nx*ny*nz];
    my=new double[nx*ny*nz];
    mz=new double[nx*ny*nz];
    en=new double[nx*ny*nz];
    bx=new double[nx*ny*nz];
    by=new double[nx*ny*nz];
    bz=new double[nx*ny*nz];
  } else{
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] ro;
    delete[] mx;
    delete[] my;
    delete[] mz;
    delete[] en;
    delete[] bx;
    delete[] by;
    delete[] bz;
    nflg=0;
  }
}
