if (n_elements(dir) eq 0) then begin
   dir=''
   read,'Set data directory: ',dir
endif

offsets=fix(reform(file_reads(dir+'/offsets.dat')))
xoff=offsets[0]
yoff=offsets[1]
zoff=offsets[2]

params=reform(file_reads(dir+'/params.dat'))
gamma=params[0]
;; beta=params[1]
;; lambda=params[2]

t=reform(file_reads(dir+'/t.dat'))
torg=t
nt=n_elements(t)

;; Read and merge X-array
fil=findfile(dir+'/x_*.dat')
mnp_x=n_elements(fil)
for i=0l,mnp_x-1 do begin
   val=file_reads(fil[i])
   nxtmp=n_elements(val)
   if (i eq 0) then x=val[xoff:nxtmp-1-xoff] else x=[x,val[xoff:nxtmp-1-xoff]]
endfor
nx=n_elements(x)

;; Read and merge Y-array
fil=findfile(dir+'/y_*.dat')
mnp_y=n_elements(fil)
for j=0l,mnp_y-1 do begin
   val=file_reads(fil[j])
   nytmp=n_elements(val)
   if (j eq 0) then y=val[yoff:nytmp-1-yoff] else y=[y,val[yoff:nytmp-1-yoff]]
endfor
ny=n_elements(y)

;; Read and merge Z-array
fil=findfile(dir+'/z_*.dat')
mnp_z=n_elements(fil)
for k=0l,mnp_z-1 do begin
   val=file_reads(fil[k])
   nztmp=n_elements(val)
   if (k eq 0) then z=val[zoff:nztmp-1-zoff] else z=[z,val[zoff:nztmp-1-zoff]]
endfor
nz=n_elements(z)

dx=x[1]-x[0]
dy=y[1]-y[0]
dz=z[1]-z[0]
idx=1./dx
idy=1./dy
idz=1./dz
dt=t[1]-t[0]

message,'Independent varaibles are loaded.',/conti
message,'To load MHD data, call READ_MHD.',/conti

end
