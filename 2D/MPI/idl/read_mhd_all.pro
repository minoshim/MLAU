;; Read MHD data in whole period

ro=fltarr(nx,ny,nt)
mx=fltarr(nx,ny,nt)
my=fltarr(nx,ny,nt)
mz=fltarr(nx,ny,nt)
en=fltarr(nx,ny,nt)
bx=fltarr(nx,ny,nt)
by=fltarr(nx,ny,nt)
bz=fltarr(nx,ny,nt)

for sst=0,nt-1 do begin

   print,"Read MHD data at "+strcompress(string(sst),/remove)+' step'
   ftmp=dir+'/outdat_'+string(sst,format='(i05)')+'*.dat'
   fil=findfile(ftmp)
   if fil[0] ne '' then begin
      for j=0l,mnp_y-1 do begin
         for i=0l,mnp_x-1 do begin
            ss=mnp_x*j+i
            val=fltarr(nxtmp,nytmp,8)
            binary_read,val,file=fil[ss]
            is=i*nx/mnp_x
            ie=(i+1)*nx/mnp_x-1
            js=j*ny/mnp_y
            je=(j+1)*ny/mnp_y-1
            ro[is:ie,js:je,sst]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,0]
            mx[is:ie,js:je,sst]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,1]
            my[is:ie,js:je,sst]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,2]
            mz[is:ie,js:je,sst]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,3]
            en[is:ie,js:je,sst]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,4]
            bx[is:ie,js:je,sst]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,5]
            by[is:ie,js:je,sst]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,6]
            bz[is:ie,js:je,sst]=val[xoff:nxtmp-1-xoff,yoff:nytmp-1-yoff,7]
         endfor
      endfor
   endif

endfor

ux=mx/ro
uy=my/ro
uz=mz/ro

;; CT B and center B
;; 2nd order for simplicity
bx_ct=bx
by_ct=by
for i=0,nx-2 do bx[i,*,*]=0.5*(bx_ct[i,*,*]+bx_ct[i+1,*,*])
for j=0,ny-2 do by[*,j,*]=0.5*(by_ct[*,j,*]+by_ct[*,j+1,*])
b2=bx^2+by^2+bz^2

pr=float((gamma-1)*(en-0.5*(ro*(ux*ux+uy*uy+uz*uz)+b2)))
te=2*pr/ro

end
