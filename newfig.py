import numpy
import matplotlib.pyplot as plt
import os




def figs(filepath, file, filoutpath, savename, fontsize=12, color='black', sfactor=1, xzoom=[-1,-1], yzoom=[-1,-1], log=0, pulsed=0):
    data = numpy.genfromtxt(filepath+file+".g2.run/g2", delimiter = ',')
        
    os.chdir(filepath+file+".g2.run/")
    count = 0 
    val = 0
    binneddata = []
    xbinneddata = []
    for row in range(len(data[:,0])):
        if data[row,0] < data[row,1]:
            count = count + 1
            val = val + data[row,4]
            if count % sfactor == 0:
                binneddata.append(val)
                xbinneddata.append(data[row,3]/1000)
                val = 0

    edit = "-0then1SMOOTH" +str(sfactor)+ "-"

    fig = plt.figure()
    
    ax = fig.add_subplot(1, 1, 1)

    ax.plot(xbinneddata,binneddata,color=color,linewidth=0.5)
    #ax.scatter(xfdata,ff,color=color, s = 0.1)#,linewidth=0.5)
    if fontsize > 12:
        ax.xaxis.set_tick_params(size = 7, width=2)
        ax.yaxis.set_tick_params(size = 7, width=2)
    
    fig.tight_layout()
    fig.savefig(savename + ".png")
    ax.set_xlim((-100,100))
    fig.savefig(savename + edit + "zoom.png")

    fdata = []
    xfdata = []
    for row in range(len(data[:,0])):
        if data[row,0] < data[row,1]:
            fdata.append(data[row,4])
            xfdata.append(data[row,3]/1000)

    edit = "-0then1-"

    fig = plt.figure()
    
    ax = fig.add_subplot(1, 1, 1)

    ax.scatter(xfdata,fdata,color=color, s = 0.1)#,linewidth=0.5)
    #ax.scatter(xfdata,ff,color=color, s = 0.1)#,linewidth=0.5)
    if fontsize > 12:
        ax.xaxis.set_tick_params(size = 7, width=2)
        ax.yaxis.set_tick_params(size = 7, width=2)
    
    fig.tight_layout()
    fig.savefig(savename + ".png")
    ax.set_xlim((-100,100))
    fig.savefig(savename + edit + "zoom.png")

    fig = plt.figure()
    
    ax = fig.add_subplot(1, 1, 1)

    ax.plot(xfdata,fdata,color=color,linewidth=0.5)
    #ax.plot(xfdata,ff,color=color,linewidth=0.5)
    if fontsize > 12:
        ax.xaxis.set_tick_params(size = 7, width=2)
        ax.yaxis.set_tick_params(size = 7, width=2)
    fig.tight_layout()
    fig.savefig(savename + edit + "line.png")
    ax.set_xlim((-100,100))
    fig.savefig(savename + edit + "linezoom.png")


    