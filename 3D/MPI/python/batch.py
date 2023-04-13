# Python3 script to load and draw MHD-3D(MPI-merged) data
# Packages Numpy and Matplotlib are required.

# Call the script in command line:
# > python batch.py
# Call the script in Python3 interactive mode:
# >>> exec(open("batch.py").read())

import numpy as np
import matplotlib.pyplot as plt
from python import plt2d

#Read independent variables and parameters
while True:
    direc=input("Input data directory (Ctrl-D to exit): ")+"/"
    try:
        x=np.loadtxt(direc+"merge_x.dat",dtype=float)
        y=np.loadtxt(direc+"merge_y.dat",dtype=float)
        z=np.loadtxt(direc+"merge_z.dat",dtype=float)
        t=np.atleast_1d(np.loadtxt(direc+"t.dat",dtype=float))
        offs=np.loadtxt(direc+"offsets.dat",dtype=int)
        para=np.atleast_1d(np.loadtxt(direc+"params.dat",dtype=float))
        break
    except:
        print("Error during file load.")

xoff=offs[0]
yoff=offs[1]
zoff=offs[2]
gam=para[0]

#Number of elements
nx=np.size(x)
ny=np.size(y)
nz=np.size(z)
nt=np.size(t)
nd=8 #Number of dependent variables in MHD-2D

#Read MHD data @ particular time
sst=-1
while ((sst < 0) or (sst >= nt)):
    sst=int(input(f"Specity time period (0-{nt-1}): "))
    
data=np.fromfile(direc+f"merge_outdat_{sst:05d}.dat",dtype=np.float32).reshape((nd,nz,ny,nx))

dx=x[1]-x[0]
dy=y[1]-y[0]
dz=z[1]-z[0]

#Primitive variables
ro=data[0,:,:,:]
vx=data[1,:,:,:]/ro
vy=data[2,:,:,:]/ro
vz=data[3,:,:,:]/ro
en=data[4,:,:,:]
bx=data[5,:,:,:]                  # @ CT grid (i-1/2,j,k)
by=data[6,:,:,:]                  # @ CT grid (i,j-1/2,k)
bz=data[7,:,:,:]                  # @ CT grid (i,j,k-1/2)
pr=(gam-1)*(en-0.5*(ro*(vx**2+vy**2+vz**2)+(bx**2+by**2+bz**2)))
