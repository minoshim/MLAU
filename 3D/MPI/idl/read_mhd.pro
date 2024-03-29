;; Read MHD data

junk=reform(file_reads(dir+'/t.dat'))
nt=n_elements(junk)

sst=0
read,'Set time step (0 - '+strcompress(string(nt-1),/remove)+') : ',sst
print,"Read MHD data at "+strcompress(string(sst),/remove)+' step'

ro=fltarr(nx,ny,nz)
mx=fltarr(nx,ny,nz)
my=fltarr(nx,ny,nz)
mz=fltarr(nx,ny,nz)
en=fltarr(nx,ny,nz)
bx=fltarr(nx,ny,nz)
by=fltarr(nx,ny,nz)
bz=fltarr(nx,ny,nz)

ftmp=dir+'/outdat_'+string(sst,format='(i05)')+'*.dat'
fil=findfile(ftmp)

if fil[0] ne '' then begin
   for k=0l,mnp_z-1 do begin
      for j=0l,mnp_y-1 do begin
         for i=0l,mnp_x-1 do begin
            ss=mnp_x*(mnp_y*k+j)+i
            val=fltarr(nxtmp,nytmp,nztmp,8)
            binary_read,val,file=fil[ss]
            is=i*nx/mnp_x
            ie=(i+1)*nx/mnp_x-1
            js=j*ny/mnp_y
            je=(j+1)*ny/mnp_y-1
            ks=k*nz/mnp_z
            ke=(k+1)*nz/mnp_z-1
            ro[is:ie,js:je,ks:ke]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,zoff:nztmp-1-zoff,0]
            mx[is:ie,js:je,ks:ke]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,zoff:nztmp-1-zoff,1]
            my[is:ie,js:je,ks:ke]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,zoff:nztmp-1-zoff,2]
            mz[is:ie,js:je,ks:ke]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,zoff:nztmp-1-zoff,3]
            en[is:ie,js:je,ks:ke]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,zoff:nztmp-1-zoff,4]
            bx[is:ie,js:je,ks:ke]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,zoff:nztmp-1-zoff,5]
            by[is:ie,js:je,ks:ke]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,zoff:nztmp-1-zoff,6]
            bz[is:ie,js:je,ks:ke]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,zoff:nztmp-1-zoff,7]
         endfor
      endfor
   endfor
endif

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

pr=float((gamma-1)*(en-0.5*(ro*(ux*ux+uy*uy+uz*uz)+b2)))
te=2*pr/ro

end
