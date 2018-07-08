import os
import subprocess
import math
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
import photon_correlation as pc

import makefig2 as m
import newfig as nm 
import newmakefigpulse as npf

import write2 as w
import analyze2 as a

def isbf(file):
    rval = 'r'
    for i in range(len(file)-1):
        if file[i] == 'b' and file[i+1] == 'f':
            rval = 'b'
            break

    return rval

def simulate(filepath, filedir, fullfilename, write = 1, analyze = 1, makefig = 1, diffuse = 1, antibunch = 1,
             pulsed = 0, endsigcts = 500000, numlines = 10**7, maxlines = 10**8, endtime = 10**12,
             temp = 298, concentration = 2*10**-8, dabsXsec = 3.6*10**-10, labsXsec = 3.6*10**(-10),
             k_demission = 1000000, k_sem = 10000, emwavelength = 815, r = 10,
             eta = 8.9 * 10**(-13), n = 1.3, k_tem = 1, k_fiss = 1, k_trans = 200000,
             reprate = 1, wavelength = 532, laserpwr = 0.5, pulselength = 80, foclen = 310000,
             NA = 1.4, darkcounts = 1, sensitivity = 0.1, nligands = 1, deadtime = 70000, afterpulse = 0, timeres = 1, order = 2,
             mode = "t2", gnpwr = 20, numbins = 4096, pulsebins = 99, channels = 3, seq = 0, mL1 = 0,
             picyzoom = 100, timestep = 200, probfiss = 1, anni = 0):
    
    os.chdir("/mnt/c/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/DotTransferSimUpdated2/")
    suffix = ".txt"
    if not os.path.isdir(filepath):
        os.mkdir(filepath)

    if not os.path.isdir(filepath + "RawData/"):
         os.mkdir(filepath + "RawData/")
         os.mkdir(filepath + "Figures/")

    if not os.path.isdir(filepath + "RawData/"+ filedir):
        os.mkdir(filepath + "RawData/"+ filedir)
        os.mkdir(filepath + "Figures/"+ filedir)
        os.mkdir(filepath + "RawData/"+ filedir+fullfilename)
        os.mkdir(filepath + "Figures/"+ filedir+fullfilename)

    elif not os.path.isdir(filepath + "RawData/"+ filedir+fullfilename):
        os.mkdir(filepath + "RawData/"+ filedir+fullfilename)
        os.mkdir(filepath + "Figures/"+ filedir+fullfilename)

    os.chdir(filepath + "RawData/"+ filedir+fullfilename)

    print(os.getcwd())
    if write == 1:
        w.write(filepath, filedir, fullfilename, antibunch, diffuse, pulsed, endsigcts, numlines, maxlines, endtime,
                    temp, concentration, dabsXsec, labsXsec,k_demission, 
                    k_fiss, k_trans, k_sem, k_tem, emwavelength,  r,
                    eta, n, reprate,wavelength, laserpwr, pulselength, foclen,
                    NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, timeres, timestep, channels, seq, mL1, probfiss, anni)

    if analyze == 1:
        a.analyze(filepath, filedir, fullfilename, numlines, order, mode, gnpwr, numbins, pulsebins, channels, makefig, m.makeafig, pulsed, picyzoom, reprate, deadtime = deadtime)


    elif makefig == 1: #only really used if it crashed part way through analysis
        #nm.figs(filepath, file, filoutpath, savename, fontsize, color, sfactor, xzoom, yzoom, log, pulsed)
        file = fullfilename
        fileoutname = file + "gnpwr" + str(gnpwr)
        if pulsed == 0:
            fileoutname = file + "gnpwr" + str(gnpwr)
            c = isbf(fullfilename)
            filename = "fig"
    
            m.makeafig("g2", filename, [-1,-1], [-1,-1], 0, pulsed, filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
            filename = "20mszoom"
            m.makeafig("g2", filename, [-20000000,20000000], [-1,-1], 0, pulsed,filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
            filename = "200uszoom"
            m.makeafig("g2", filename, [-200000,200000], [-1,-1], 0, pulsed,filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
            filename = "log"
            m.makeafig("g2", filename,[-1,-1],[-1,-1], 1, pulsed,filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)

        else:
            dfilepath = filepath + "RawData/" + filedir + file + "/" + file + "gnpwr" + str(gnpwr) + ".g2.run/"
            sfilepath = filepath + "Figures/" + filedir + file + "/"
            npf.makepulsedfig(dfilepath, "g2", sfilepath, fileoutname, deadtime = deadtime, reprate = reprate, timespace = 1000)
            npf.makepulsedfig(dfilepath, "g2", sfilepath, fileoutname+"250nszoom", deadtime = deadtime, reprate = reprate, timespace = 1000, xzoom = 100)
            npf.makepulsedfig(dfilepath, "g2", sfilepath, fileoutname+"sq", deadtime = deadtime, reprate = reprate, timespace = 1000, yzoom = [0,80], figsize = [6,4])
            npf.makepulsedfig(dfilepath, "g2", sfilepath, fileoutname+"250nszoomsq", deadtime = deadtime, reprate = reprate, timespace = 1000, xzoom = 100,yzoom = [0,80], figsize = [6,4])
'''
        print(filepath+"/"+file)
        fileout = file + "-PIC"
        filename = "PIC"
        m.makeafig(fileout, filename,[-1,-1], [-1,-1], 1, pulsed,filepath+"/"+file, color = c)
        filename = filename + "xzoom"
        m.makeafig(fileout, filename, [0,10**6], picyzoom, 1, pulsed,filepath+"/"+file, color = c)
        filename = "PICtight"
        m.makeafig(fileout,filename,[0.0001,10], picyzoom, 1, pulsed,filepath+"/"+file, color = c)
        filename = "PICmed"
        m.makeafig(fileout,filename,[0,10**10], picyzoom, 1, pulsed,filepath+"/"+file, color = c)'''




