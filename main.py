import sim as sim

#control parameters - 1 means on, 0 means off
write = 1
analyze = 1
makefig = 1

diffuse = 1
antibunch = 1
pulsed = 1

#General saving folder parameters
filepath = "/mnt/c/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/DotTransferSim/"
filedir = "May17LongfixedDotem/"
filenames = ["LowLP"]

#file length - goes to the max of these two
endsigcts = 10000000
numlines = 10**8
maxlines = 10**9
endtime = 300*10**12 #ps #10**15 = 16.6 min

#Sample parameters
temp = 298 #K
k_demission = 1400000 #emission lifetime in ps
k_sem = 10000 #short lived singlets
k_tem = 100000000 #long lived triplets
k_fiss = 1 #made up from marks paper on pentacene fission
k_trans = 200000 #based on Rao paper 2018
emwavelength = 815
r = 20 #nm - hydrodynamic radius of particles

eta = 8.9 * 10**(-13) # kg/nm s - dynamic viscosity of solvent (water)
n = 1.3 # index of refraction - sample in water

concentration = 2*10**(-8)
nligands = 100

dabsXsec = 7.180616773321853*10**(-9)# per emitter numabs'd = phperpulse*absXsec*numEms - this is reasonable based on absXsec for CdSe is 550000*r^3/cm (from Bawendi paper Ruvim sent me)
labsXsec = dabsXsec*10

#Laser parameters
reprate = 0.05 #MHz

wavelength = 532 #nm
beamdiam = 5000000 #5 mm in nm
laserpwr = 0.5 # mW (0.52 mWinto back of objective 23 #mW)
#laserpwr = [0.05,0.1,0.2,0.3,0.5,0.75,1,10]#[0.001,0.01,0.025,

pulselength = 80 #ps - not used


#Objective parameters
foclen = 310000 #310 microns in nm (working distance + coverslip thickness)
NA = 1.4

#Detector parameters
darkcounts = 100 #s^-1
sensitivity = 0.1
deadtime = 70000 #70 ns in ps
afterpulse = 0.0001 #percent of time a photon is emitted a deadtime after one is detected

#Correlation parameters
order = 2 #g2
mode = "t2"
gnpwr = 23
numbins = 8192 # 2**23/8192 = 1024 => 1 bin is 1 ns (ish)
pulsebins = (4)-1#should always be an odd number
channels = 2

#Miscellaneous simulation parameters
picyzoom = [-1,-1]
timestep = 1000000 #average number of photons per "round" of calculations

seq = 1
i = 0
f = filenames[0]
count = 0
for laserpwr in [0.00005,0.005,0.000005]:
    for concentration in [2*10**(-8),2*10**(-6),2*10**(-10)]:
        i = i + 1
        for sensitivity in [0.1,0.5,0.01]:
            for k_demission in [1400000,15]:
                for nligands in [1,100,10000]:
                    for diffuse in range(2):
                        count = count + 1
                        if count > 3:
                            f = filenames[0]
                            if i == 2:
                                f = f + "conc"
                            elif i == 3:
                                f = f + "dilute"
                            if k_demission == 15:
                                f = f + "CdSe"
                            else:
                                f = f + "PbS"
                            f = f + str(nligands) + "ligs"
                            if diffuse == 1:
                                f = f + "DIFF"
                            sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                                pulsed, endsigcts, numlines, maxlines, endtime,temp, concentration,
                                                dabsXsec, labsXsec, k_demission, k_sem,
                                                emwavelength, r,eta, n, k_tem, k_fiss, k_trans,
                                                reprate, wavelength, laserpwr, pulselength, foclen,
                                                NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, order,
                                                mode, gnpwr, numbins, pulsebins, channels, seq)


'''
for mL1 in range(2):
    for seq in range(2):
        for nl in range(2):
            nligands = 1
            if nl == 0:
                nligands = 100
            for diffuse in range(2):
                f = filenames[0] + "-seq-" + str(seq) + "-" + str(nligands)+"lig"

                if mL1 == 1:
                    f = f + "-mL1"
                if diffuse == 1:
                    f = f + "DIFF"


                sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                        pulsed, numlines, maxlines, endtime,temp, concentration,
                                        dabsXsec, labsXsec, k_demission, k_sem,
                                        emwavelength, r,eta, n, k_tem, k_fiss, k_trans,
                                        reprate, wavelength, laserpwr, pulselength, foclen,
                                        NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, order,
                                        mode, gnpwr, numbins, pulsebins, channels, seq, mL1)


f = filenames[0]
sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                    pulsed, numlines, maxlines, endtime,temp, concentration,
                                    dabsXsec, labsXsec, k_demission, k_sem,
                                    emwavelength, r,eta, n, k_tem, k_fiss, k_trans,
                                    reprate, wavelength, laserpwr, pulselength, foclen,
                                    NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, order,
                                    mode, gnpwr, numbins, pulsebins, channels, seq)

for seq in range(1):
    nligands = 100
    f = filenames[0] + "-seq-" + str(seq) + "-" + str(nligands)+"lig"
    sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                pulsed, numlines, maxlines, endtime,temp, concentration,
                                dabsXsec, labsXsec, k_demission, k_sem,
                                emwavelength, r,eta, n, k_tem, k_fiss, k_trans,
                                reprate, wavelength, laserpwr, pulselength, foclen,
                                NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, order,
                                mode, gnpwr, numbins, pulsebins, channels, seq)

    nligands = 1
    f = filenames[0] + "-seq-" + str(seq) + "mL-1lig"
    sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                    pulsed, numlines, maxlines, endtime,temp, concentration,
                                    dabsXsec, labsXsec, k_demission, k_sem,
                                    emwavelength, r,eta, n, k_tem, k_fiss, k_trans,
                                    reprate, wavelength, laserpwr, pulselength, foclen,
                                    NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, order,
                                    mode, gnpwr, numbins, pulsebins, channels, seq, mL1 = 1)

sim.simulate(filepath, filedir, "2018-03-27-shortish-seq-1-100ligshortg2", 0, analyze, makefig, diffuse, antibunch,
                                    pulsed, numlines, maxlines, endtime,temp, concentration,
                                    dabsXsec, labsXsec, k_demission, k_sem,
                                    emwavelength, r,eta, n, k_tem, k_fiss, k_trans,
                                    reprate, wavelength, laserpwr, pulselength, foclen,
                                    NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, order,
                                    mode, gnpwr, numbins, pulsebins, channels, seq, mL1 = 1)

for i in range(3):
    sensitivity = 0.01/8
    laserpwr = laserpwri*10**(i-1)
    for j in range(3):
        sensitivity = sensitivity*8
        for k in range(3):
            k_demission = k_demissioni*10**(k-1)

            filename = "sens-" + str(int(100*sensitivity)) + "-lp-"+str(int(1000*laserpwr)) + "uW-tauDotEm-"+str(int(k_demission/1000)) +"ns"

            '''
