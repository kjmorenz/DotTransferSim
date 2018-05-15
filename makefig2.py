import os
import subprocess
import photon_correlation as pc 

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