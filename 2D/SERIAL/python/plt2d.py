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
