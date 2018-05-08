import os
import subprocess
import math
import matplotlib.pyplot as plt
import photon_correlation as pc

import makefig2 as m
import newfig as nm 

import write as w
import analyze2 as a

def isbf(file):
    rval = 'r'
    for i in range(len(file)-1):
        if file[i] == 'b' and file[i+1] == 'f':
            rval = 'b'
            break

    return rval

def simulate(filepath, filedir, fullfilename, write = 1, analyze = 1, makefig = 1, diffuse = 1, antibunch = 1,
             pulsed = 0, numlines = 10**7, maxlines = 10**8, endtime = 10**12,
             temp = 298, concentration = 2*10**-8, dabsXsec = 3.6*10**-10, labsXsec = 3.6*10**(-10),
             k_demission = 1000000, k_sem = 10000, emwavelength = 815, r = 10,
             eta = 8.9 * 10**(-13), n = 1.3, k_tem = 1, k_fiss = 1, k_trans = 200000,
             reprate = 1, wavelength = 532, laserpwr = 0.5, pulselength = 80, foclen = 310000,
             NA = 1.4, darkcounts = 1, sensitivity = 0.1, nligands = 1, deadtime = 70000, afterpulse = 0, order = 2,
             mode = "t2", gnpwr = 20, numbins = 4096, pulsebins = 99, channels = 3, seq = 0, mL1 = 0,
             picyzoom = 100, timestep = 200):
    
    
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
        w.write(filepath, filedir, fullfilename, antibunch, diffuse, pulsed, numlines, maxlines, endtime,
                    temp, concentration, dabsXsec, labsXsec,k_demission, 
                    k_fiss, k_trans, k_sem, k_tem, emwavelength,  r,
                    eta, n, reprate,wavelength, laserpwr, pulselength, foclen,
                    NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, timestep, channels, seq, mL1)

    if analyze == 1:
        a.analyze(filepath, filedir, fullfilename, numlines, order, mode, gnpwr, numbins, pulsebins, channels, makefig, m.makeafig, pulsed, picyzoom, reprate)


    elif makefig == 1: #only really used if it crashed part way through analysis
        #nm.figs(filepath, file, filoutpath, savename, fontsize, color, sfactor, xzoom, yzoom, log, pulsed)

        file = fullfilename
        fileoutname = file + "gnpwr15"
        c = isbf(fullfilename)
        filename = "fig"
  
        m.makeafig("g2", filename, [-1,-1], [-1,-1], 0, pulsed, filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
        filename = "20mszoom"
        m.makeafig("g2", filename, [-20000000,20000000], [-1,-1], 0, pulsed,filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
        filename = "200uszoom"
        m.makeafig("g2", filename, [-200000,200000], [-1,-1], 0, pulsed,filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
        filename = "log"
        m.makeafig("g2", filename,[-1,-1],[-1,-1], 1, pulsed,filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
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




