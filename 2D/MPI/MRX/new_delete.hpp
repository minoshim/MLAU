void new_delete()
{
  static int nflg=0;
  if (nflg == 0){
    nflg=1;
    x=new double[nx];
    y=new double[ny];
    ro=new double[nx*ny];
    mx=new double[nx*ny];
    my=new double[nx*ny];
    mz=new double[nx*ny];
    en=new double[nx*ny];
    bx=new double[nx*ny];
    by=new double[nx*ny];
    bz=new double[nx*ny];
  } else{
    delete[] x;
    delete[] y;
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
