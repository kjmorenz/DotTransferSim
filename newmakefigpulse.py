from io import BytesIO
import numpy
import matplotlib.pyplot as plt
import matplotlib

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

def tickfunction(trep, npulses):
    ticks = []
    labels = []
    for i in range(2*npulses + 1):
        ticks.append(trep*(i-npulses))
        labels.append(str(i-npulses))
    return ticks, labels

def makepulsedfig(filepath, file, savename, reprate = 1, npulses=1, time=float('inf'), fontsize=12, figsize = [-1,-1]):
    trep = (10**6)/reprate # in ps
    data, xdata, bindata, cdata, cxdata, endbin = effload(filepath,file,100000,npulses,time,trep)

    fig = plt.figure()
    fig.patch.set_facecolor('black')
    ax = fig.add_subplot(1, 1, 1)
    ax2 = ax.twiny()

    
    matplotlib.rcParams.update({'font.size': fontsize})
    matplotlib.rc('xtick', labelsize = fontsize)

    color = 'white'

    ax.plot(xdata,data, color)
    ax.plot(cxdata,cdata, color = 'red')
    ticks, labels = tickfunction(trep, npulses)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ticks)
    ax2.set_xticklabels(labels)
    ax.set_axis_bgcolor('none')
    ax.spines['top'].set_color('white')
    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')

    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')

    ax2.set_axis_bgcolor('none')
    ax2.spines['top'].set_color('white')
    ax2.spines['bottom'].set_color('white')
    ax2.spines['left'].set_color('white')
    ax2.spines['right'].set_color('white')

    ax2.xaxis.label.set_color('white')
    ax2.yaxis.label.set_color('white')

    ax.set_ylabel("$g^{(2)}$")
    ax.set_xlabel("time relative to pulse (ns)")
    ax2.set_xlabel("pulses elapsed")

    if fontsize > 12:
        ax.xaxis.set_tick_params(size = 7, width=2)
        ax.yaxis.set_tick_params(size = 7, width=2)

    if not figsize == [-1,-1]:
        fig.set_size_inches(figsize[0], figsize[1])

    ax.tick_params(axis = 'x', colors = 'white')
    ax.tick_params(axis = 'y', colors = 'white')


    fig.tight_layout()
    fig.savefig(filepath + savename,facecolor=fig.get_facecolor(), edgecolor = 'none')
    plt.close()