if (n_elements(dir) eq 0) then begin
   dir=''
   read,'Set data directory: ',dir
endif

offsets=fix(reform(file_reads(dir+'/offsets.dat')))
xoff=offsets[0]
yoff=offsets[1]
zoff=offsets[2]

x=reform(file_reads(dir+'/x.dat'))
y=reform(file_reads(dir+'/y.dat'))
z=reform(file_reads(dir+'/z.dat'))
t=reform(file_reads(dir+'/t.dat'))
gamma=(reform(file_reads(dir+'/params.dat')))[0]
nx=n_elements(x)
ny=n_elements(y)
nz=n_elements(z)
nt=n_elements(t)
dx=x[1]-x[0]
dy=y[1]-y[0]
dz=z[1]-z[0]

x=x[xoff:nx-1-xoff]
y=y[yoff:ny-1-yoff]
z=z[zoff:nz-1-zoff]

message,'Independent varaibles are loaded.',/conti
message,'To load MHD data, call READ_MHD.',/conti
end

