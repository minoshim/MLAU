void setbc()
{
  // Set boundary condition
  bc3d(ro,nx,ny,nz,XOFF,YOFF,ZOFF,0,0,0,0,0,0);
  bc3d(mx,nx,ny,nz,XOFF,YOFF,ZOFF,0,0,0,0,0,0);
  bc3d(my,nx,ny,nz,XOFF,YOFF,ZOFF,0,0,0,0,0,0);
  bc3d(mz,nx,ny,nz,XOFF,YOFF,ZOFF,0,0,0,0,0,0);
  bc3d(en,nx,ny,nz,XOFF,YOFF,ZOFF,0,0,0,0,0,0);
  bc3d(bx,nx,ny,nz,XOFF,YOFF,ZOFF,1,0,0,0,0,0);
  bc3d(by,nx,ny,nz,XOFF,YOFF,ZOFF,0,0,1,0,0,0);
  bc3d(bz,nx,ny,nz,XOFF,YOFF,ZOFF,0,0,0,0,1,0);
}
