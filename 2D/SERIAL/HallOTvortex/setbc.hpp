void setbc()
{
  // Set boundary condition
  bc2d(ro,nx,ny,XOFF,YOFF,0,0,0,0);
  bc2d(mx,nx,ny,XOFF,YOFF,0,0,0,0);
  bc2d(my,nx,ny,XOFF,YOFF,0,0,0,0);
  bc2d(mz,nx,ny,XOFF,YOFF,0,0,0,0);
  bc2d(en,nx,ny,XOFF,YOFF,0,0,0,0);
  bc2d(bx,nx,ny,XOFF,YOFF,1,0,0,0);
  bc2d(by,nx,ny,XOFF,YOFF,0,0,1,0);
  bc2d(bz,nx,ny,XOFF,YOFF,0,0,0,0);
}
