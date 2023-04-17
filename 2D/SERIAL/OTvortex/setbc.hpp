void setbc()
{
  // Set boundary condition
  double* p[8]={ro,mx,my,mz,en,bx,by,bz};
  for (int n=0;n<8;n++) bc2d(p[n],nx,ny,xoff,yoff,stxs[n],dnxs[n],stys[n],dnys[n]);
}
