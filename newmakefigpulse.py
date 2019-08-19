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
def findafterglow(deadtime, cxdata):
    diffpos = float("inf")
    diffneg = float("inf")
    indexpos = 0
    indexneg = 0
    for i in range(len(cxdata)):
        if numpy.abs(cxdata[i]-deadtime) < diffpos:
            secondindpos = indexpos
            diffpos = numpy.abs(cxdata[i]-deadtime)
            indexpos = i
        if numpy.abs(cxdata[i] + deadtime) < diffneg:
            secondindneg = indexneg
            diffneg = numpy.abs(cxdata[i]+deadtime)
            indexneg = i
    return indexpos, secondindpos, indexneg, secondindneg
    #print(cxdata[indexpos])
    #print(cxdata[indexneg])

def effload(filepath,file,printevery,npulses,time,trep, g2norm =1):
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
            if thisline.size > 4:
                if val == 1 and numpy.absolute(thisline[4]) < time:
                    endbin = numpy.absolute(thisline[4])
                    
                if thisline[0] < thisline[1] and max(numpy.absolute(thisline[2]), numpy.absolute(thisline[3])) <= npulses + 0.5 and numpy.absolute(thisline[4]) < time:
                    
                    bindata.append(thisline[2])
                    if thisline[2] == -0.5:
                        cdata.append(thisline[-1]/g2norm)
                        cxdata.append(thisline[4]+(thisline[2]+0.5)*trep)
                    else:
                        data.append(thisline[-1]/g2norm)
                        xdata.append(thisline[4]+(thisline[2]+0.5)*trep)
    return data, xdata, bindata, cdata, cxdata, endbin

def tickfunction(trep, npulses, divby):
    ticks = []
    labels = []
    for i in range(2*npulses + 1):
        ticks.append(trep*(i-npulses)/divby)
        labels.append(str(i-npulses))
    return ticks, labels

def multimakepulsedfig(filepath, file, sfilepath, savenames, reprate = 1, npulses=1, time=float('inf'), fontsize=12, figsizes = [[-1,-1]], linewidth = 2, 
                    grayafterglow = 1, deadtime = 70000, timespace = -1, xzooms = [-1], yzooms =[[-1,-1]], bbg = 1, g2norm = 1):
    
    
    if bbg == 1:
        spcolor = 'white'
        agcolor = 'darkslategray'
        bbgname = 'bbg'
    elif bbg == 0:
        spcolor = 'black'
        agcolor = "gray"
        bbgname = 'wbg'
    cpcolor = 'red'

    trep = (10**6)/reprate # in ps
    data, xdata, bindata, cdata, cxdata, endbin = effload(filepath,file,10000,npulses,time,trep, g2norm)
    if grayafterglow ==1:
        inp, sinp, inn, sinn = findafterglow(deadtime, cxdata)

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

    if grayafterglow ==1:
        if len(cxdata) > inp + 2 and len(cxdata) > inn + 2:
            agxpdata = [cxdata[inp-1], cxdata[inp], cxdata[inp+1]]
            agxndata = [cxdata[inn-1], cxdata[inn], cxdata[inn+1]]
            agpdata = [cdata[inp-1], cdata[inp], cdata[inp+1]]
            agndata = [cdata[inn-1], cdata[inn], cdata[inn+1]]
            del cxdata[inp-1], cxdata[inp-1], cxdata[inp-1], cxdata[inn-1], cxdata[inn-1], cxdata[inn-1]
            del cdata[inp-1], cdata[inp-1], cdata[inp-1], cdata[inn-1], cdata[inn-1], cdata[inn-1]
        else: grayafterglow = 0

    for i in range(len(savenames)):
        savename = savenames[i]
        figsize = figsizes[i]
        xzoom = xzooms[i]
        yzoom = yzooms[i]
      
        fig = plt.figure()
        if bbg == 1:
            fig.patch.set_facecolor('black')
        ax = fig.add_subplot(1, 1, 1)
        ax2 = ax.twiny()

        
        matplotlib.rcParams.update({'font.size': fontsize})
        matplotlib.rc('xtick', labelsize = fontsize)

        ax.plot(xdata[:len(cxdata)],data[:len(cxdata)], color = spcolor, linewidth = linewidth)
        ax.plot(xdata[len(cxdata)+7:],data[len(cxdata)+7:], color = spcolor, linewidth = linewidth)
        
        if grayafterglow == 1:
            ax.plot(agxpdata,agpdata, color = agcolor, linewidth = linewidth)
            ax.plot(agxndata,agndata, color = agcolor, linewidth = linewidth)
            
        
        ax.plot(cxdata,cdata, color = cpcolor, linewidth = linewidth)
        

        
        ticks, labels = tickfunction(trep, npulses, divby)
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(ticks)
        ax2.set_xticklabels(labels)
        if bbg == 1:
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

        plt.ylabel("$g^{(2)}($" + u'\u03C4' + "$)$")
        plt.xlabel("$" + u'\u03C4' + "$ (" + timetag + ")" )
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
        
        if bbg == 1:
            ax.tick_params(axis = 'x', colors = 'white')
            ax.tick_params(axis = 'y', colors = 'white')
            ax2.tick_params(axis = 'x', colors = 'white')
            ax2.tick_params(axis = 'y', colors = 'white')


        fig.tight_layout()
        fig.savefig(sfilepath + savename + bbgname + ".png",facecolor=fig.get_facecolor(), edgecolor = 'none')
        
        plt.close()

        #symmetric

        fig = plt.figure()
        if bbg == 1:
            fig.patch.set_facecolor('black')
        ax = fig.add_subplot(1, 1, 1)
        ax2 = ax.twiny()

        
        matplotlib.rcParams.update({'font.size': fontsize})
        matplotlib.rc('xtick', labelsize = fontsize)

        ax.plot(xdata[:len(cxdata)],data[:len(cxdata)], color = spcolor, linewidth = linewidth)
        ax.plot(xdata[len(cxdata)+7:],data[len(cxdata)+7:], color = spcolor, linewidth = linewidth)
        for i in range(len(xdata)):
            xdata[i] = -xdata[i]
        ax.plot(xdata[:len(cxdata)],data[:len(cxdata)], color = spcolor, linewidth = linewidth)
        ax.plot(xdata[len(cxdata)+7:],data[len(cxdata)+7:], color = spcolor, linewidth = linewidth)
        if grayafterglow == 1:
            ax.plot(agxpdata,agpdata, color = agcolor, linewidth = linewidth)
            ax.plot(agxndata,agndata, color = agcolor, linewidth = linewidth)
            for i in range(len(agxpdata)):
                agxpdata[i] = - agxpdata[i]
                agxndata[i] = - agxndata[i]
            ax.plot(agxpdata,agpdata, color = agcolor, linewidth = linewidth)
            ax.plot(agxndata,agndata, color = agcolor, linewidth = linewidth)
        ax.plot(cxdata,cdata, color = cpcolor, linewidth = linewidth)
        for i in range(len(cxdata)):
            cxdata[i] = -cxdata[i]
        ax.plot(cxdata,cdata, color = cpcolor, linewidth = linewidth)
        

        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(ticks)
        ax2.set_xticklabels(labels)
        if bbg ==1:
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

        plt.ylabel("$g^{(2)}($" + u'\u03C4' + "$)$")
        plt.xlabel("$" + u'\u03C4' + "$ (" + timetag + ")" )
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
        
        if bbg == 1:
            ax.tick_params(axis = 'x', colors = 'white')
            ax.tick_params(axis = 'y', colors = 'white')
            ax2.tick_params(axis = 'x', colors = 'white')
            ax2.tick_params(axis = 'y', colors = 'white')


        fig.tight_layout()
        fig.savefig(sfilepath + savename + bbgname + "SYMM.png",facecolor=fig.get_facecolor(), edgecolor = 'none')
        
        plt.close()

        if xzoom != -1:
            ax.set_xlim((-xzoom, xzoom))
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(ticks)
            ax2.set_xticklabels(labels)
            if yzoom != [-1,-1]:
                ax.set_ylim((yzoom[0],yzoom[1]))
                sname = savename + "SYMM-xzoom" +str(xzoom)+"yzoom" + str(yzoom[0])+"to"+str(yzoom[1]) +  bbgname + ".png"
            else:
                sname = savename + "SYMM-zoom" +str(xzoom)+ bbgname + ".png"
            if not figsize == [-1,-1]:
                fig.set_size_inches(figsize[0], figsize[1])
            else:
                fig.set_size_inches(6, 4)

            fig.savefig(sfilepath + sname, facecolor=fig.get_facecolor(), edgecolor = 'none')
            plt.close()
        elif yzoom !=[-1,-1]:
            ax.set_ylim((yzoom[0],yzoom[1]))
            fig.savefig(sfilepath + savename + "SYMMyzoom" + str(yzoom[0])+"to"+str(yzoom[1]) +  bbgname + ".png",facecolor=fig.get_facecolor(), edgecolor = 'none')
            plt.close()
        

        #replot overlap

        #symmetric
        linewidth = linewidth/2
        fig = plt.figure()
        if bbg == 1:
            fig.patch.set_facecolor('black')
        ax = fig.add_subplot(1, 1, 1)
        
        matplotlib.rcParams.update({'font.size': fontsize})
        matplotlib.rc('xtick', labelsize = fontsize)

        ax.plot(cxdata,data[:len(cxdata)], color = spcolor, linewidth = linewidth)
        for i in range(len(cxdata)):
            cxdata[i] = -cxdata[i]
        ax.plot(cxdata,data[:len(cxdata)], color = spcolor, linewidth = linewidth)
        if grayafterglow == 1:
            ax.plot(agxpdata,agpdata, color = agcolor, linewidth = linewidth)
            ax.plot(agxndata,agndata, color = agcolor, linewidth = linewidth)
            for i in range(len(agxpdata)):
                agxpdata[i] = - agxpdata[i]
                agxndata[i] = - agxndata[i]
            ax.plot(agxpdata,agpdata, color = agcolor, linewidth = linewidth)
            ax.plot(agxndata,agndata, color = agcolor, linewidth = linewidth)
        ax.plot(cxdata,cdata, color = cpcolor, linewidth = linewidth)
        for i in range(len(cxdata)):
            cxdata[i] = -cxdata[i]
        ax.plot(cxdata,cdata, color = cpcolor, linewidth = linewidth)
        

        if bbg == 1:
            ax.set_facecolor('black')
            ax.spines['top'].set_color('white')
            ax.spines['bottom'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.spines['right'].set_color('white')

            ax.xaxis.label.set_color('white')
            ax.yaxis.label.set_color('white')

        

        plt.ylabel("$g^{(2)}($" + u'\u03C4' + "$)$")
        plt.xlabel("$" + u'\u03C4' + "$ (" + timetag + ")" )


        if fontsize > 12:
            ax.xaxis.set_tick_params(size = 7, width=2)
            ax.yaxis.set_tick_params(size = 7, width=2)

        if not figsize == [-1,-1]:
            fig.set_size_inches(figsize[0], figsize[1])
        else:
            fig.set_size_inches(6, 4)
        
        if bbg == 1:
            ax.tick_params(axis = 'x', colors = 'white')
            ax.tick_params(axis = 'y', colors = 'white')


        fig.tight_layout()
        fig.savefig(sfilepath + savename + "OVERLAPSYMM" + bbgname + ".png",facecolor=fig.get_facecolor(), edgecolor = 'none')
        
        plt.close()

        if xzoom != -1:
            ax.set_xlim((-xzoom, xzoom))
            if yzoom != [-1,-1]:
                ax.set_ylim((yzoom[0],yzoom[1]))
                savename = savename + "OVERLAPSYMM-xzoom" +str(xzoom)+"yzoom" + str(yzoom[0])+"to"+str(yzoom[1]) +  bbgname + ".png"
            else:
                savename = savename + "OVERLAPSYMM-zoom" +str(xzoom)+ bbgname +".png"
            if not figsize == [-1,-1]:
                fig.set_size_inches(figsize[0], figsize[1])
            else:
                fig.set_size_inches(5, 4)
            fig.savefig(sfilepath + savename, facecolor=fig.get_facecolor(), edgecolor = 'none')
            plt.close()
        elif yzoom !=[-1,-1]:
            ax.set_ylim((yzoom[0],yzoom[1]))
            fig.savefig(sfilepath + savename + "OVERLAPSYMMyzoom" + str(yzoom[0])+"to"+str(yzoom[1]) +  bbgname + ".png",facecolor=fig.get_facecolor(), edgecolor = 'none')
            plt.close()
        
def makepulsedfig(filepath, file, sfilepath, savename, reprate = 1, npulses=1, time=float('inf'), fontsize=12, figsize = [-1,-1], linewidth = 0.5, 
                    grayafterglow = 1, deadtime = 70000, timespace = -1, xzoom = -1, yzoom =[-1,-1], bbg = 1):
    if bbg == 1:
        agcolor = "darkslategray"
        spcolor = 'white'
        cpcolor = 'red'
    elif bbg == 0:
        agcolor = "gray"
        spcolor = 'black'
        cpcolor = 'red'

    trep = (10**6)/reprate # in ps
    data, xdata, bindata, cdata, cxdata, endbin = effload(filepath,file,10000,npulses,time,trep)
    if grayafterglow ==1:
        inp, sinp, inn, sinn = findafterglow(deadtime, cxdata)
    
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
    if grayafterglow ==1:
        if len(cxdata) > inp + 2 and len(cxdata) > inn + 2:
            agxpdata = [cxdata[inp-1], cxdata[inp], cxdata[inp+1]]
            agxndata = [cxdata[inn-1], cxdata[inn], cxdata[inn+1]]
            agpdata = [cdata[inp-1], cdata[inp], cdata[inp+1]]
            agndata = [cdata[inn-1], cdata[inn], cdata[inn+1]]
            del cxdata[inp-1], cxdata[inp-1], cxdata[inp-1], cxdata[inn-1], cxdata[inn-1], cxdata[inn-1]
            del cdata[inp-1], cdata[inp-1], cdata[inp-1], cdata[inn-1], cdata[inn-1], cdata[inn-1]
        else: grayafterglow = 0

    fig = plt.figure()
    if bbg == 1:
        fig.patch.set_facecolor('black')

    ax = fig.add_subplot(1, 1, 1)
    ax2 = ax.twiny()

    
    matplotlib.rcParams.update({'font.size': fontsize})
    matplotlib.rc('xtick', labelsize = fontsize)
    
    ax.plot(xdata[:len(cxdata)],data[:len(cxdata)], color = spcolor, linewidth = linewidth)
    ax.plot(xdata[len(cxdata)+7:],data[len(cxdata)+7:], color = spcolor, linewidth = linewidth)

    if grayafterglow == 1:
        ax.plot(agxpdata,agpdata, color = agcolor, linewidth = linewidth)
        ax.plot(agxndata,agndata, color = agcolor, linewidth = linewidth)
        
    
    ax.plot(cxdata,cdata, color = cpcolor, linewidth = linewidth)
    

    
    ticks, labels = tickfunction(trep, npulses, divby)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ticks)
    ax2.set_xticklabels(labels)
    if bbg == 1:
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

    plt.ylabel("$g^{(2)}($" + u'\u03C4' + "$)$")
    plt.xlabel("$" + u'\u03C4' + "$ (" + timetag + ")" )
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
    if bbg == 1:
        ax.tick_params(axis = 'x', colors = 'white')
        ax.tick_params(axis = 'y', colors = 'white')
        ax2.tick_params(axis = 'x', colors = 'white')
        ax2.tick_params(axis = 'y', colors = 'white')


    fig.tight_layout()
    fig.savefig(sfilepath + savename + ".png",facecolor=fig.get_facecolor(), edgecolor = 'none')
    
    plt.close()

    #symmetric

    fig = plt.figure()
    if bbg == 1:
        fig.patch.set_facecolor('black')
    ax = fig.add_subplot(1, 1, 1)
    ax2 = ax.twiny()

    
    matplotlib.rcParams.update({'font.size': fontsize})
    matplotlib.rc('xtick', labelsize = fontsize)

    ax.plot(xdata[:len(cxdata)],data[:len(cxdata)], color = spcolor, linewidth = linewidth)
    ax.plot(xdata[len(cxdata)+7:],data[len(cxdata)+7:], color = spcolor, linewidth = linewidth)
    for i in range(len(xdata)):
        xdata[i] = -xdata[i]
    ax.plot(xdata[:len(cxdata)],data[:len(cxdata)], color = spcolor, linewidth = linewidth)
    ax.plot(xdata[len(cxdata)+7:],data[len(cxdata)+7:], color = spcolor, linewidth = linewidth)
    if grayafterglow == 1:
        ax.plot(agxpdata,agpdata, color = agcolor, linewidth = linewidth)
        ax.plot(agxndata,agndata, color = agcolor, linewidth = linewidth)
        for i in range(len(agxpdata)):
            agxpdata[i] = - agxpdata[i]
            agxndata[i] = - agxndata[i]
        ax.plot(agxpdata,agpdata, color = agcolor, linewidth = linewidth)
        ax.plot(agxndata,agndata, color = agcolor, linewidth = linewidth)
    ax.plot(cxdata,cdata, color = cpcolor, linewidth = linewidth)
    for i in range(len(cxdata)):
        cxdata[i] = -cxdata[i]
    ax.plot(cxdata,cdata, color = cpcolor, linewidth = linewidth)
    

    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ticks)
    ax2.set_xticklabels(labels)
    if bbg == 1:
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

    plt.ylabel("$g^{(2)}($" + u'\u03C4' + "$)$")
    plt.xlabel("$" + u'\u03C4' + "$ ("+timetag+")",size=fontsize)
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
    
    if bbg == 1:
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
            sname = savename + "SYMM-xzoom" +str(xzoom)+"yzoom" + str(yzoom[0])+"to"+str(yzoom[1]) +  ".png"
        else:
            sname = savename + "SYMM-zoom" +str(xzoom)+".png"
        if not figsize == [-1,-1]:
            fig.set_size_inches(figsize[0], figsize[1])
        else:
            fig.set_size_inches(6, 4)
        fig.savefig(sfilepath + sname, facecolor=fig.get_facecolor(), edgecolor = 'none')
        plt.close()
    elif yzoom !=[-1,-1]:
        ax.set_ylim((yzoom[0],yzoom[1]))
        fig.savefig(sfilepath + savename + "SYMMyzoom" + str(yzoom[0])+"to"+str(yzoom[1]) +  ".png",facecolor=fig.get_facecolor(), edgecolor = 'none')
        plt.close()
    

    #replot overlap

    #symmetric
    linewidth = linewidth/2
    fig = plt.figure()
    if bbg == 1:
        fig.patch.set_facecolor('black')
    ax = fig.add_subplot(1, 1, 1)
    
    matplotlib.rcParams.update({'font.size': fontsize})
    matplotlib.rc('xtick', labelsize = fontsize)

    ax.plot(cxdata,data[:len(cxdata)], color = spcolor, linewidth = linewidth)
    for i in range(len(cxdata)):
        cxdata[i] = -cxdata[i]
    ax.plot(cxdata,data[:len(cxdata)], color = spcolor, linewidth = linewidth)
    if grayafterglow == 1:
        ax.plot(agxpdata,agpdata, color = agcolor, linewidth = linewidth)
        ax.plot(agxndata,agndata, color = agcolor, linewidth = linewidth)
        for i in range(len(agxpdata)):
            agxpdata[i] = - agxpdata[i]
            agxndata[i] = - agxndata[i]
        ax.plot(agxpdata,agpdata, color = agcolor, linewidth = linewidth)
        ax.plot(agxndata,agndata, color = agcolor, linewidth = linewidth)
    ax.plot(cxdata,cdata, color = cpcolor, linewidth = linewidth)
    for i in range(len(cxdata)):
        cxdata[i] = -cxdata[i]
    ax.plot(cxdata,cdata, color = cpcolor, linewidth = linewidth)
    

    if bbg == 1:
        ax.set_facecolor('black')
        ax.spines['top'].set_color('white')
        ax.spines['bottom'].set_color('white')
        ax.spines['left'].set_color('white')
        ax.spines['right'].set_color('white')

        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')

    

    plt.ylabel("$g^{(2)}($" + u'\u03C4' + "$)$")
    plt.xlabel("$" + u'\u03C4' + "$ ("+timetag+")",size=fontsize)


    if fontsize > 12:
        ax.xaxis.set_tick_params(size = 7, width=2)
        ax.yaxis.set_tick_params(size = 7, width=2)

    if not figsize == [-1,-1]:
        fig.set_size_inches(figsize[0], figsize[1])
    else:
        fig.set_size_inches(6, 4)
    
    if bbg == 1:
        ax.tick_params(axis = 'x', colors = 'white')
        ax.tick_params(axis = 'y', colors = 'white')


    fig.tight_layout()
    fig.savefig(sfilepath + savename + "OVERLAPSYMM.png",facecolor=fig.get_facecolor(), edgecolor = 'none')
    
    plt.close()

    if xzoom != -1:
        ax.set_xlim((-xzoom, xzoom))
        if yzoom != [-1,-1]:
            ax.set_ylim((yzoom[0],yzoom[1]))
            savename = savename + "OVERLAPSYMM-xzoom" +str(xzoom)+"yzoom" + str(yzoom[0])+"to"+str(yzoom[1]) +  ".png"
        else:
            savename = savename + "OVERLAPSYMM-zoom" +str(xzoom)+".png"
        if not figsize == [-1,-1]:
            fig.set_size_inches(figsize[0], figsize[1])
        else:
            fig.set_size_inches(5, 4)
        fig.savefig(sfilepath + savename, facecolor=fig.get_facecolor(), edgecolor = 'none')
        plt.close()
    elif yzoom !=[-1,-1]:
        ax.set_ylim((yzoom[0],yzoom[1]))
        fig.savefig(sfilepath + savename + "OVERLAPSYMMyzoom" + str(yzoom[0])+"to"+str(yzoom[1]) +  ".png",facecolor=fig.get_facecolor(), edgecolor = 'none')
        plt.close()
'''    
filepath = "/mnt/c/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/DotTransferSim/"
filedir = "May18LongfixedDotem/"
file = "LowLPno10kligsPbS1ligs"
gnpwr = 23
fileoutname = file
dfilepath = filepath + "RawData/" + filedir + file + "/" + file + "gnpwr" + str(gnpwr) + ".g2.run/"
sfilepath = filepath + "Figures/" + filedir + file + "/"
deadtime = 70000
reprate = 0.0312

multimakepulsedfig(dfilepath, "g2", sfilepath, [fileoutname+"---", fileoutname+"2500nszoom",fileoutname+"sq",fileoutname+"2500nszoomsq"], 
                                deadtime = deadtime, reprate = reprate, timespace = 1000000, g2norm = 100,
                                xzooms = [0,2500,0,2500], yzooms = [[-1,-1],[-1,-1],[0,20],[0,20]], 
                                figsizes = [[-1,-1], [-1,-1],[6,4],[6,4]], bbg = 1, fontsize = 20)'''