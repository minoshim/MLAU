import matplotlib.pyplot as plt

def image(x,y,val,save=0,title="",cmap="jet",xlabel="x",ylabel="y",figsize=(6.4,4.8),filename="result.eps",equal=0,show=0,aspect="auto"):
    fig=plt.figure(figsize=figsize)

    if (equal != 0):
        aspect="equal"
    plt.axes().set_aspect(aspect)

    # plt.pcolormesh(x,y,val,cmap=cmap,shading='auto')
    plt.imshow(val,cmap=cmap,
               extent=[x.min(),x.max(),y.min(),y.max()],
               interpolation='nearest',origin='lower',aspect=aspect) # Fast

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.colorbar()
    if (save):
        fig.savefig(filename)
    if (show != 0):
        plt.show(block=False) # console non-blocked
    return fig

def imagex3(x,y,z,val,save=0,title="",cmap="jet"):
    fig=plt.figure(figsize=(9,3))
    nx=len(x)
    ny=len(y)
    nz=len(z)

    plt.subplot(1,3,1)
    tmp=0.5*(val[nz//2-1,:,:]+val[nz//2,:,:])
    plt.pcolormesh(x,y,tmp,cmap=cmap,shading='auto')
    plt.xlabel("x")
    plt.ylabel("y")

    plt.subplot(1,3,2)
    tmp=0.5*(val[:,:,nx//2-1]+val[:,:,nx//2])
    plt.pcolormesh(y,z,tmp,cmap=cmap,shading='auto')
    plt.xlabel("y")
    plt.ylabel("z")
    plt.title(title)
    
    plt.subplot(1,3,3)
    tmp=0.5*(val[:,ny//2-1,:]+val[:,ny//2,:])
    plt.pcolormesh(x,z,tmp,cmap=cmap,shading='auto')
    plt.xlabel("x")
    plt.ylabel("z")

    fig.tight_layout()
    plt.show(block=False) # console non-blocked
    return fig
