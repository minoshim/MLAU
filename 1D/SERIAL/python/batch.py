# Python3 script to load and draw MHD-1D data
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
        x=np.loadtxt(direc+"x.dat",dtype=np.float64)
        t=np.loadtxt(direc+"t.dat",dtype=np.float64)
        xoff=int(np.loadtxt(direc+"xoff.dat",dtype=int))
        para=np.loadtxt(direc+"params.dat",dtype=np.float64)
        break
    except:
        print("Error during file load.")
    
bx0=para[0]
gam=para[1]

#Number of elements
nx=np.size(x)
nt=np.size(t)
nd=7 #Number of dependent variables in MHD-1D

#Read MHD data @ particular time
sst=-1
while ((sst < 0) or (sst >= nt)):
    sst=int(input(f"Specity time period (0-{nt-1}): "))
    
data=np.fromfile(direc+f"outdat_{sst:05d}.dat",dtype=np.float64).reshape((nd,nx))

#Slice to remove ghost cells
x=x[xoff:nx-xoff]
data2=data[:,xoff:nx-xoff]

#Primitive variables
ro=data2[0,:]
vx=data2[1,:]/ro
vy=data2[2,:]/ro
vz=data2[3,:]/ro
bx=np.array([bx0]*(nx-2*xoff))
by=data2[5,:]
bz=data2[6,:]
pr=(gam-1)*(data2[4,:]-0.5*(ro*(vx**2+vy**2+vz**2)+(bx**2+by**2+bz**2)))
data2=np.array([ro,vx,vy,vz,pr,bx,by,bz])

#Plot
ylabels=("Density","Vx","Vy","Vz","Pressure","Bx","By","Bz")
fig=plt.figure(figsize=(10,4))
fig.suptitle(f"t={t[sst]:.3f}")
for i in range(0,8):
    plt.subplot(2,4,i+1)
    plt.plot(x,data2[i,:])
    plt.xlabel("x")
    plt.ylabel(ylabels[i])

fig.tight_layout()
fig.savefig("result.eps")
# plt.show() #console blocked
plt.show(block=False) # console non-blocked
