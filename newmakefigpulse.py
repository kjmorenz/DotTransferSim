from io import BytesIO
import numpy
import matplotlib
matplotlib.use('Agg')#FOR LINUX, OTHERWISE AUTOBACKEND IS XWINDOWS???

'''otherwise gives error:
File "main.py", line 106, in <module>
    mode, gnpwr, numbins, pulsebins, channels, seq)
  File "/mnt/c/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/DotTransferSimUpdated2/sim.py", line 73, in simulate
    m.makeafig("g2", filename, [-1,-1], [-1,-1], 0, pulsed, filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
  File "/mnt/c/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/DotTransferSimUpdated2/makefig2.py", line 32, in makeafig
    fig = g2.make_figure(log, xzoom = xzoom, yzoom = yzoom, fontsize = fontsize, normalize = normalize, scale = scale)#if log = 0 then it's not logscale, 1 is logscale
  File "/mnt/c/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/DotTransferSimUpdated2/photon_correlation/G2.py", line 129, in make_figure
    fig = plt.figure()
  File "/home/karen/.local/lib/python3.5/site-packages/matplotlib/pyplot.py", line 548, in figure
    **kwargs)
  File "/home/karen/.local/lib/python3.5/site-packages/matplotlib/backend_bases.py", line 161, in new_figure_manager
    return cls.new_figure_manager_given_figure(num, fig)
  File "/home/karen/.local/lib/python3.5/site-packages/matplotlib/backends/_backend_tk.py", line 1044, in new_figure_manager_given_figure
    window = Tk.Tk(className="matplotlib")
  File "/usr/lib/python3.5/tkinter/__init__.py", line 1871, in __init__
    self.tk = _tkinter.create(screenName, baseName, className, interactive, wantobjects, useTk, sync, use)
_tkinter.TclError: no display name and no $DISPLAY environment variable'''
import matplotlib.pyplot as plt

def effload(filepath,file,printevery,npulses,time,trep):
    val = 0
    data = []
    xdata = []
    bindata = []
    cdata = []
    cxdata = []
    endbin = time
    with open(filepath+file, 'r') as fileobj:
        for line in fileobj:
            val = val + 1
            if printevery != 0 and val%printevery == 0:
                print(str(val) + " lines")
            thisline = numpy.genfromtxt(BytesIO(line.encode('utf-8')), delimiter = ",")
            if val == 1 and numpy.absolute(thisline[4]) < time:
                endbin = numpy.absolute(thisline[4])
                
            if thisline[0] < thisline[1] and max(numpy.absolute(thisline[2]), numpy.absolute(thisline[3])) <= npulses + 0.5 and numpy.absolute(thisline[4]) < time:
                data.append(thisline[-1])
                xdata.append(thisline[4]+(thisline[2]+0.5)*trep)
                bindata.append(thisline[2])
                if thisline[2] == -0.5:
                    cdata.append(thisline[-1])
                    cxdata.append(thisline[4]+(thisline[2]+0.5)*trep)
    return data, xdata, bindata, cdata, cxdata, endbin

def tickfunction(trep, npulses, divby):
    ticks = []
    labels = []
    for i in range(2*npulses + 1):
        ticks.append(trep*(i-npulses)/divby)
        labels.append(str(i-npulses))
    return ticks, labels

def makepulsedfig(filepath, file, sfilepath, savename, reprate = 1, npulses=1, time=float('inf'), fontsize=12, figsize = [-1,-1], linewidth = 0.5, timespace = -1, xzoom = -1, yzoom =[-1,-1]):
    trep = (10**6)/reprate # in ps
    data, xdata, bindata, cdata, cxdata, endbin = effload(filepath,file,100000,npulses,time,trep)
    divby = 1
    if not timespace == -1:
        if 1000 <= timespace <= 100000:
            timetag = "ns"
            divby = 1000
        elif 100000 <= timespace <= 100000000:
            timetag = "us"
            divby = 1000000
        elif 100000000 <= timespace <= 100000000000:
            timetag = "ms"
            divby = 1000000000
        else:
            timetag = "s"
            divby = 10**12
        for i in range(len(xdata)):
            xdata[i] = xdata[i]/divby
        for i in range(len(cxdata)):
            cxdata[i] = cxdata[i]/divby


    fig = plt.figure()
    fig.patch.set_facecolor('black')
    ax = fig.add_subplot(1, 1, 1)
    ax2 = ax.twiny()

    
    matplotlib.rcParams.update({'font.size': fontsize})
    matplotlib.rc('xtick', labelsize = fontsize)

    color = 'white'

    ax.plot(xdata,data, color = color, linewidth = linewidth)
    ax.plot(cxdata,cdata, color = 'red', linewidth = linewidth)
    
    ticks, labels = tickfunction(trep, npulses, divby)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ticks)
    ax2.set_xticklabels(labels)
    ax.set_facecolor('black')
    ax.spines['top'].set_color('white')
    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')

    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')

    ax2.set_facecolor('black')
    ax2.spines['top'].set_color('white')
    ax2.spines['bottom'].set_color('white')
    ax2.spines['left'].set_color('white')
    ax2.spines['right'].set_color('white')

    ax2.xaxis.label.set_color('white')
    ax2.yaxis.label.set_color('white')

    ax.set_ylabel("$g^{(2)}$")
    ax.set_xlabel("time relative to pulse (" + timetag + ")" )
    ax2.set_xlabel("pulses elapsed")

    if fontsize > 12:
        ax.xaxis.set_tick_params(size = 7, width=2)
        ax.yaxis.set_tick_params(size = 7, width=2)
        ax2.xaxis.set_tick_params(size = 7, width=2)
        ax2.yaxis.set_tick_params(size = 7, width=2)

    if not figsize == [-1,-1]:
        fig.set_size_inches(figsize[0], figsize[1])
    else:
        fig.set_size_inches(3*(2*npulses + 1), 4)
    ax.tick_params(axis = 'x', colors = 'white')
    ax.tick_params(axis = 'y', colors = 'white')
    ax2.tick_params(axis = 'x', colors = 'white')
    ax2.tick_params(axis = 'y', colors = 'white')


    fig.tight_layout()
    fig.savefig(sfilepath + savename + ".png",facecolor=fig.get_facecolor(), edgecolor = 'none')
    
    plt.close()

    #symmetric

    fig = plt.figure()
    fig.patch.set_facecolor('black')
    ax = fig.add_subplot(1, 1, 1)
    ax2 = ax.twiny()

    
    matplotlib.rcParams.update({'font.size': fontsize})
    matplotlib.rc('xtick', labelsize = fontsize)

    color = 'white'

    ax.plot(xdata,data, color = color, linewidth = linewidth)
    ax.plot(cxdata,cdata, color = 'red', linewidth = linewidth)
    for i in range(len(xdata)):
        xdata[i] = -xdata[i]
    for i in range(len(cxdata)):
        cxdata[i] = -cxdata[i]
    
    ax.plot(xdata,data, color = color, linewidth = linewidth)
    ax.plot(cxdata,cdata, color = 'red', linewidth = linewidth)

    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ticks)
    ax2.set_xticklabels(labels)
    ax.set_facecolor('black')
    ax.spines['top'].set_color('white')
    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')

    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')

    ax2.set_facecolor('black')
    ax2.spines['top'].set_color('white')
    ax2.spines['bottom'].set_color('white')
    ax2.spines['left'].set_color('white')
    ax2.spines['right'].set_color('white')

    ax2.xaxis.label.set_color('white')
    ax2.yaxis.label.set_color('white')

    ax.set_ylabel("$g^{(2)}$")
    ax.set_xlabel("time relative to pulse (" + timetag + ")" )
    ax2.set_xlabel("pulses elapsed")

    if fontsize > 12:
        ax.xaxis.set_tick_params(size = 7, width=2)
        ax.yaxis.set_tick_params(size = 7, width=2)
        ax2.xaxis.set_tick_params(size = 7, width=2)
        ax2.yaxis.set_tick_params(size = 7, width=2)

    if not figsize == [-1,-1]:
        fig.set_size_inches(figsize[0], figsize[1])
    else:
        fig.set_size_inches(3*(2*npulses + 1), 4)
    ax.tick_params(axis = 'x', colors = 'white')
    ax.tick_params(axis = 'y', colors = 'white')
    ax2.tick_params(axis = 'x', colors = 'white')
    ax2.tick_params(axis = 'y', colors = 'white')


    fig.tight_layout()
    fig.savefig(sfilepath + savename + "SYMM.png",facecolor=fig.get_facecolor(), edgecolor = 'none')
    
    plt.close()

    if xzoom != -1:
        ax.set_xlim((-xzoom, xzoom))
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(ticks)
        ax2.set_xticklabels(labels)
        if yzoom != [-1,-1]:
            ax.set_ylim((yzoom[0],yzoom[1]))
            savename = savename + "SYMM-xzoom" +str(xzoom)+"yzoom" + str(yzoom[0])+"to"+str(yzoom[1]) +  ".png"
        else:
            savename = savename + "SYMM-zoom" +str(xzoom)+".png"
        if not figsize == [-1,-1]:
            fig.set_size_inches(figsize[0], figsize[1])
        else:
            fig.set_size_inches(6, 4)
        fig.savefig(sfilepath + savename, facecolor=fig.get_facecolor(), edgecolor = 'none')
        plt.close()
    elif yzoom !=[-1,-1]:
        ax.set_ylim((yzoom[0],yzoom[1]))
        fig.savefig(sfilepath + savename + "SYMMyzoom" + str(yzoom[0])+"to"+str(yzoom[1]) +  ".png",facecolor=fig.get_facecolor(), edgecolor = 'none')
        plt.close()
'''filepath = '/mnt/c/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/DotTransferSim/RawData/May15-test/fixedDeadtime70psto70nsPbS100ligs/fixedDeadtime70psto70nsPbS100ligsgnpwr23.g2.run/'
file = "g2"
savename = '/mnt/c/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/DotTransferSim/Figures/May15-test/fixedDeadtime70psto70nsPbS100ligs/practicenewpulsedfig'
reprate = 0.05
time = 2**23
makepulsedfig(filepath, file, filepath, savename, reprate = reprate, timespace = 1000, xzoom = 100)
'''