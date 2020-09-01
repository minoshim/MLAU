# Python3 script to load and draw MHD-2D(Serial) data
# Packages Numpy and Matplotlib are required.

# Call the script in command line:
# > python batch.py
# Call the script in Python3 interactive mode:
# >>> exec(open("batch.py").read())

import numpy as np
import matplotlib.pyplot as plt

#Read independent variables and parameters
while True:
    direc=input("Input data directory (Ctrl-D to exit): ")+"/"
    try:
        x=np.loadtxt(direc+"x.dat",dtype=np.float)
        y=np.loadtxt(direc+"y.dat",dtype=np.float)
        t=np.loadtxt(direc+"t.dat",dtype=np.float)
        offs=np.loadtxt(direc+"offsets.dat",dtype=np.int)
        para=float(np.loadtxt(direc+"params.dat",dtype=np.float))
        break
    except:
        print("Error during file load.")

xoff=offs[0]
yoff=offs[1]
gam=para

#Number of elements
nx=np.size(x)
ny=np.size(y)
nt=np.size(t)
nd=8 #Number of dependent variables in MHD-2D

#Read MHD data @ particular time
sst=-1
while ((sst < 0) or (sst >= nt)):
    sst=int(input(f"Specity time period (0-{nt-1}): "))
    
data=np.fromfile(direc+f"outdat_{sst:05d}.dat",dtype=np.float64).reshape((nd,nx,ny))

#Slice to remove ghost cells
x=x[xoff:nx-xoff]
y=y[yoff:ny-yoff]
data2=data[:,xoff:nx-xoff,yoff:ny-yoff]

#Primitive variables
ro=data2[0,:,:]
vx=data2[1,:,:]/ro
vy=data2[2,:,:]/ro
vz=data2[3,:,:]/ro
bx=data2[5,:,:]
by=data2[6,:,:]
bz=data2[7,:,:]
pr=(gam-1)*(data2[4,:,:]-0.5*(ro*(vx**2+vy**2+vz**2)+(bx**2+by**2+bz**2)))
data2=np.array([ro,vx,vy,vz,pr,bx,by,bz])

#Plot
val=pr/ro
fig=plt.figure()
plt.pcolormesh(x,y,val,cmap="jet")
plt.xlabel("x")
plt.ylabel("y")
plt.title(f"t={t[sst]:.2f}")
plt.colorbar()
fig.savefig("result.eps")
plt.show(block=False) # console non-blocked
