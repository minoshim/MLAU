if (n_elements(dir) eq 0) then begin
   dir=''
   read,'Set data directory: ',dir
endif

xoff=fix((reform(file_reads(dir+'/xoff.dat')))[0])
bx=(reform(file_reads(dir+'/params.dat')))[0]
gamma=(reform(file_reads(dir+'/params.dat')))[1]

x=reform(file_reads(dir+'/x.dat'))
t=reform(file_reads(dir+'/t.dat'))
nx=n_elements(x)
nt=n_elements(t)
dx=x[1]-x[0]
dt=t[1]-t[0]

files=findfile(dir+'/outdat_*.dat')
ro=dblarr(nx,nt)
mx=dblarr(nx,nt)
my=dblarr(nx,nt)
mz=dblarr(nx,nt)
en=dblarr(nx,nt)
by=dblarr(nx,nt)
bz=dblarr(nx,nt)
for i=0,n_elements(files)-1 do begin
   tmp=dblarr(nx,7)
   binary_read,tmp,file=files[i]
   ro[*,i]=tmp[*,0]
   mx[*,i]=tmp[*,1]
   my[*,i]=tmp[*,2]
   mz[*,i]=tmp[*,3]
   en[*,i]=tmp[*,4]
   by[*,i]=tmp[*,5]
   bz[*,i]=tmp[*,6]
endfor
bx=replicate(bx,nx,nt)

x=x[xoff:nx-1-xoff]
ro=ro[xoff:nx-1-xoff,*]
mx=mx[xoff:nx-1-xoff,*]
my=my[xoff:nx-1-xoff,*]
mz=mz[xoff:nx-1-xoff,*]
en=en[xoff:nx-1-xoff,*]
bx=bx[xoff:nx-1-xoff,*]
by=by[xoff:nx-1-xoff,*]
bz=bz[xoff:nx-1-xoff,*]
nx=n_elements(x)

vx=mx/ro
vy=my/ro
vz=mz/ro
pr=(gamma-1)*(en-0.5*(ro*(vx*vx+vy*vy+vz*vz)+(bx*bx+by*by+bz*bz)))

end
