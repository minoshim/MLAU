if (n_elements(dir) eq 0) then begin
   dir=''
   read,'Set data directory: ',dir
endif

offsets=fix(reform(file_reads(dir+'/offsets.dat')))
xoff=offsets[0]
yoff=offsets[1]

x=reform(file_reads(dir+'/x.dat'))
y=reform(file_reads(dir+'/y.dat'))
t=reform(file_reads(dir+'/t.dat'))
gamma=(reform(file_reads(dir+'/params.dat')))[0]
nx=n_elements(x)
ny=n_elements(y)
nt=n_elements(t)
dx=x[1]-x[0]
dy=y[1]-y[0]
dt=t[1]-t[0]

files=findfile(dir+'/outdat_*.dat')
ro=dblarr(nx,ny,nt)
mx=dblarr(nx,ny,nt)
my=dblarr(nx,ny,nt)
mz=dblarr(nx,ny,nt)
en=dblarr(nx,ny,nt)
bx=dblarr(nx,ny,nt)
by=dblarr(nx,ny,nt)
bz=dblarr(nx,ny,nt)
for i=0,n_elements(files)-1 do begin
   tmp=dblarr(nx,ny,8)
   binary_read,tmp,file=files[i]
   ro[*,*,i]=tmp[*,*,0]
   mx[*,*,i]=tmp[*,*,1]
   my[*,*,i]=tmp[*,*,2]
   mz[*,*,i]=tmp[*,*,3]
   en[*,*,i]=tmp[*,*,4]
   bx[*,*,i]=tmp[*,*,5]
   by[*,*,i]=tmp[*,*,6]
   bz[*,*,i]=tmp[*,*,7]
endfor

x=x[xoff:nx-1-xoff]
y=y[yoff:ny-1-yoff]
ro=ro[xoff:nx-1-xoff,yoff:ny-1-yoff,*]
mx=mx[xoff:nx-1-xoff,yoff:ny-1-yoff,*]
my=my[xoff:nx-1-xoff,yoff:ny-1-yoff,*]
mz=mz[xoff:nx-1-xoff,yoff:ny-1-yoff,*]
en=en[xoff:nx-1-xoff,yoff:ny-1-yoff,*]
bx=bx[xoff:nx-1-xoff,yoff:ny-1-yoff,*]
by=by[xoff:nx-1-xoff,yoff:ny-1-yoff,*]
bz=bz[xoff:nx-1-xoff,yoff:ny-1-yoff,*]
nx=n_elements(x)
ny=n_elements(y)

vx=mx/ro
vy=my/ro
vz=mz/ro
pr=(gamma-1)*(en-0.5*(ro*(vx*vx+vy*vy+vz*vz)+(bx*bx+by*by+bz*bz)))

end
