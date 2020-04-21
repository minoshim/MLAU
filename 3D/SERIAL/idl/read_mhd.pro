junk=reform(file_reads(dir+'/t.dat'))
nt=n_elements(junk)

sst=0
read,'Set time step (0 - '+strcompress(string(nt-1),/remove)+') : ',sst
print,"Read MHD data at "+strcompress(string(sst),/remove)+' step'
val=dblarr(nx,ny,nz,8)
binary_read,val,file=dir+'/outdat_'+string(sst,format='(i05)')+'.dat'
ro=val[*,*,*,0]
mx=val[*,*,*,1]
my=val[*,*,*,2]
mz=val[*,*,*,3]
en=val[*,*,*,4]
bx=val[*,*,*,5]
by=val[*,*,*,6]
bz=val[*,*,*,7]

ux=mx/ro
uy=my/ro
uz=mz/ro

;; CT B and center B
;; 2nd order for simplicity
bx_ct=bx
by_ct=by
bz_ct=bz
for i=0,nx-2 do bx[i,*,*]=0.5*(bx_ct[i,*,*]+bx_ct[i+1,*,*])
for j=0,ny-2 do by[*,j,*]=0.5*(by_ct[*,j,*]+by_ct[*,j+1,*])
for k=0,nz-2 do bz[*,*,k]=0.5*(bz_ct[*,*,k]+bz_ct[*,*,k+1])
b2=bx^2+by^2+bz^2

pr=(gamma-1)*(en-0.5*(ro*(ux*ux+uy*uy+uz*uz)+b2))

ro=ro[xoff:nx-1-xoff,yoff:ny-1-yoff,zoff:nz-1-zoff]
mx=mx[xoff:nx-1-xoff,yoff:ny-1-yoff,zoff:nz-1-zoff]
my=my[xoff:nx-1-xoff,yoff:ny-1-yoff,zoff:nz-1-zoff]
mz=mz[xoff:nx-1-xoff,yoff:ny-1-yoff,zoff:nz-1-zoff]
ux=ux[xoff:nx-1-xoff,yoff:ny-1-yoff,zoff:nz-1-zoff]
uy=uy[xoff:nx-1-xoff,yoff:ny-1-yoff,zoff:nz-1-zoff]
uz=uz[xoff:nx-1-xoff,yoff:ny-1-yoff,zoff:nz-1-zoff]
en=en[xoff:nx-1-xoff,yoff:ny-1-yoff,zoff:nz-1-zoff]
pr=pr[xoff:nx-1-xoff,yoff:ny-1-yoff,zoff:nz-1-zoff]
bx=bx[xoff:nx-1-xoff,yoff:ny-1-yoff,zoff:nz-1-zoff]
by=by[xoff:nx-1-xoff,yoff:ny-1-yoff,zoff:nz-1-zoff]
bz=bz[xoff:nx-1-xoff,yoff:ny-1-yoff,zoff:nz-1-zoff]

te=2*pr/ro

end
