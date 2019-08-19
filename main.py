import sim as sim

#Run in Ubuntu not cygwin because rust not installed on cygwin

#control parameters - 1 means on, 0 means off
write = 1
analyze = 1
makefig = 1

diffuse = 1
antibunch = 1
pulsed = 1

#General saving folder parameters
filepath = "/mnt/c/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/DotTransferSim/"
filedir = "2019-07-26-numems/"
filenames = ["Stationary fewem-low-diffusing"]
filename = filenames[0]

#file length - goes to the max of these two
endsigcts = 100000000
numlines = 10**8#5
maxlines = 10**9
endtime = 2*(10**5)*10**12 #ps #10**15 = 16.6 min

#Sample parameters
temp = 298 #K
k_demission = 1400000 #emission lifetime in ps
k_sem = 10000 #short lived singlets #1000000000000 # never singlet emit from Yb #
k_tem = 100000000 #long lived triplets #1000000000000 # pretend Yb has unity qy, never 'triplet emits' #
k_fiss = 1 #made up from marks paper on pentacene fission
k_trans = 200000# [50000, 200000, 5000000]  #based on Rao paper 2018 1#
emwavelength = 815
r = 20 #nm - hydrodynamic radius of particles

eta = 8.9 * 10**(-13) # kg/nm s - dynamic viscosity of solvent (water)
n = 1.3 # index of refraction - sample in water

concentration = 2*10**(-8)#2*10**(-8) #I think there aren't many Yb so I left this at 'low concentration'
concs = [1,2,10,3]
nligands = [1] #i.e. absorbtion centres which transfer to Yb

dabsXsec = 7.180616773321853*10**(-9)# per emitter numabs'd = phperpulse*absXsec*numEms - this is reasonable based on absXsec for CdSe is 550000*r^3/cm (from Bawendi paper Ruvim sent me)
labsXsec = dabsXsec*10

#Laser parameters
reprate = 0.05#312 #MHz

wavelength = 532 #nm 470 #
beamdiam = 5000000 #5 mm in nm
laserpwr = 5*10**(-4) # mW (0.52 mWinto back of objective 23 #mW)
laserpwr = [0.001, 0.01, 0.05,0.1,0.2,0.3,0.5,0.75,1,10]#[0.001,0.01,0.025,
laserpwr = [10**(-4), 10**(-5), 10**(-6), 10**(-7), 10**(-8), 10**(-9)]#, 0.001, 0.01, 0.05,0.1,0.2,0.3,0.5,0.75,1,10]
laserpwr = 5*10**(-5)
pulselength = 80 #ps - not used


#Objective parameters
foclen = 310000 #310 microns in nm (working distance + coverslip thickness)
NA = 1.4

#Detector parameters
darkcounts = 100 #s^-1
sensitivity = 0.01
deadtime = 10000000#70000 #70 ns in ps
deadtime = 70000
afterpulse = 0.0001 #percent of time a photon is emitted a deadtime after one is detected
timeres = 400 #Not Used

#Correlation parameters
order = 2 #g2
mode = "t2"
gnpwr = 23
numbins = 8192 # 2**23/8192 = 1024 => 1 bin is 1 ns (ish)
pulsebins = (4)-1#should always be an odd number
channels = 2

#Miscellaneous simulation parameters
picyzoom = [-1,-1]
timestep = 5000 #average number of photons per "round" of calculations



'''
def simulate(filepath, filedir, fullfilename, write = 1, analyze = 1, makefig = 1, diffuse = 1, antibunch = 1,
             pulsed = 0, endsigcts = 500000, numlines = 10**7, maxlines = 10**8, endtime = 10**12,
             temp = 298, concentration = 2*10**-8, dabsXsec = 3.6*10**-10, labsXsec = 3.6*10**(-10),
             k_demission = 1000000, k_sem = 10000, emwavelength = 815, r = 10,
             eta = 8.9 * 10**(-13), n = 1.3, k_tem = 1, k_fiss = 1, k_trans = 200000,
             reprate = 1, wavelength = 532, laserpwr = 0.5, pulselength = 80, foclen = 310000,
             NA = 1.4, darkcounts = 1, sensitivity = 0.1, nligands = 1, deadtime = 70000, afterpulse = 0, timeres = 1, order = 2,
             mode = "t2", gnpwr = 20, numbins = 4096, pulsebins = 99, channels = 3, seq = 0, mL1 = 0,
             picyzoom = 100, timestep = 200, probfiss = 1, anni = 0):
'''

for nems in concs:
    f = filename + str(nems) + "ems"
    print(f)
    sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                    pulsed, endsigcts, numlines, maxlines, endtime,temp, concentration*nems,
                    dabsXsec, labsXsec, k_demission, k_sem,
                    emwavelength, r,eta, n, k_tem, k_fiss, k_trans,
                    reprate, wavelength, laserpwr, pulselength, foclen,
                    NA, darkcounts, sensitivity, 1, deadtime, afterpulse, timeres, order,
                    mode, gnpwr, numbins, pulsebins, channels, timestep=timestep,  printall = 0)


'''
for k_trans in [k_sem]:#k_transs:
    f = filenames[0] + str(int(k_trans/1000)) + "ns-10ligs-lowlp"
    sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                    pulsed, endsigcts, numlines, maxlines, endtime,temp, concentration,
                    dabsXsec, labsXsec, k_demission, k_sem,
                    emwavelength, r,eta, nligands, k_tem, k_fiss, k_trans,
                    reprate, wavelength, laserpwr, pulselength, foclen,
                    NA, darkcounts, sensitivity, 1, deadtime, afterpulse, timeres, order,
                    mode, gnpwr, numbins, pulsebins, channels, timestep=timestep)

count = 0
for dots in ["PbS"]:
    for diffuse in range(2):
        for nligands in [10000]:
            count = count + 1
            
            f =  filenames[0] + dots + str(nligands) + "ligs"
            
            if diffuse == 1:
                f= f + "DIFF"
            #f = f + "sens1"
            if count ==1 or count ==2:
                sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                pulsed, endsigcts, numlines, maxlines, endtime,temp, concentration,
                                dabsXsec, labsXsec/nligands, k_demission, k_sem,
                                emwavelength, r,eta, n, k_tem, k_fiss, k_trans,
                                reprate, wavelength, laserpwr, pulselength, foclen,
                                NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, timeres, order,
                                mode, gnpwr, numbins, pulsebins, channels, timestep=timestep)
        




for anni in [0]:
    for diffuse in range(1,2):
        for lpf in [10, 1]:
            
            f =  filenames[0]
            if anni == 1:
                f= f + "anni"
            if diffuse == 1:
                f= f + "DIFF"
                endsigcts = 5*10**7
            else:
                endsigcts = 10**7
            if lpf == 10:
                f = f + "highlp"
            sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                            pulsed, endsigcts, numlines, maxlines, endtime,temp, concentration,
                            dabsXsec, labsXsec/nligands, k_demission, k_sem,
                            emwavelength, r,eta, n, k_tem, k_fiss, k_trans,
                            reprate, wavelength, laserpwr*lpf, pulselength, foclen,
                            NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, timeres, order,
                            mode, gnpwr, numbins, pulsebins, channels, timestep=timestep, anni=anni)
        
seq = 1
i = 0
f = filenames[0]
anni = 0
#for probfiss in [0]:
for anni in range(1):
    for laserpwr in [0.00005]:
        for k_demission in [1400000, 15000]:
            for nligands in [1,100]:
                for diffuse in range(2):
                    f = filenames[0]
                    if k_demission == 15000:
                        f = f + "CdSe"
                    else:
                        f = f + "PbS"
                    f = f + str(nligands) + "ligs"
                    if diffuse == 1:
                        f = f + "DIFF"
                    if anni == 1:
                        f = f + "anni"
                    #if probfiss == 0:
                    #    f = f + "NOFISS"
                    sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                        pulsed, endsigcts, numlines, maxlines, endtime,temp, concentration,
                                        dabsXsec, labsXsec/nligands, k_demission, k_sem,
                                        emwavelength, r,eta, n, k_tem, k_fiss, k_trans,
                                        reprate, wavelength, laserpwr, pulselength, foclen,
                                        NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, order,
                                        mode, gnpwr, numbins, pulsebins, channels, seq, anni = anni)

filedir = "May28actualoldexparam/"
filenames = [""]
f = filenames[0]
count = 0
for probfiss in [1,0]:
    for anni in range(2):
        for laserpwr in [0.00005]:
            for k_demission in [1400000, 15000]:
                for nligands in [1,100]:
                    for diffuse in range(2):
                        if nligands == 100:
                            count =1
                        if count == 1:
                            f = filenames[0]
                            if k_demission == 15000:
                                f = f + "CdSe"
                            else:
                                f = f + "PbS"
                            f = f + str(nligands) + "ligs"
                            if diffuse == 1:
                                f = f + "DIFF"
                            if anni == 1:
                                f = f + "anni"
                            if probfiss == 0:
                                f = f + "NOFISS"
                            sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                                pulsed, endsigcts, numlines, maxlines, endtime,temp, concentration,
                                                dabsXsec, 10*dabsXsec, k_demission, k_sem,
                                                emwavelength, r,eta, n, k_tem, k_fiss, k_trans,
                                                reprate, wavelength, laserpwr, pulselength, foclen,
                                                NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, order,
                                                mode, gnpwr, numbins, pulsebins, channels, seq, probfiss=probfiss, anni = anni)'''

