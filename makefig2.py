import os
import subprocess
import photon_correlation as pc 
import matplotlib.pyplot as plt

'''
Make fig takes g2 files and makes plots out of them with input parameters
g2file - the name of the file with the g2 correlation in it (usually "g2")
filename - the name to save the figure under
xzoom - a tuple with the x min and max values to plot
yzoom - a tuple with the y min and max values to plot
log - 1 if you want logscale x axis, 0 otherwise
pulsed - 1 if your data is pulsed, 0 otherwise
fontsize - size of font on plot
filepath - if you aren't sitting in the directory with the g2 file, tell the program where the g2 file is located
normalize - the value to divide all your g2 bins by
scale - a yshift to add to all your bins
color - the color of the trace on the plot
'''
def makeafig(g2file,filename, xzoom, yzoom, log, pulsed, fontsize = 12, filepath = "",filedir = "",fileoutdir = "", normalize = 1, scale = 0, color = 'r'):
    if pulsed != 1:
        #print((filepath + "RawData/" + filedir + fileoutdir + g2file))
        g2 = pc.G2_T2(filepath + "RawData/" + filedir + fileoutdir + g2file) #, int_counts=False)
        
    else:
        g2 = pc.G2_T3(filepath + "RawData/" + filedir + fileoutdir+ g2file)
        '''fig = g2.make_figure(log, xzoom = xzoom, yzoom = yzoom)#if log = 0 then it's not logscale, 1 is logscale
        fig.savefig(filename+".png")
        print(fig)
        plt.clf()'''
    
    fig = g2.make_figure(log, xzoom = xzoom, yzoom = yzoom, fontsize = fontsize, normalize = normalize, scale = scale)#if log = 0 then it's not logscale, 1 is logscale
    print(fig)
    fig.savefig(filepath + "Figures/" + filedir + filename+".png")
    print(filename)
    plt.close()

def isbf(file):
    rval = 'r'
    for i in range(len(file)-1):
        if file[i] == 'b' and file[i+1] == 'f':
            rval = 'b'
            break

    return rval