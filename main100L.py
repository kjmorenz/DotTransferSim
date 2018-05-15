import sim as sim

#control parameters - 1 means on, 0 means off
write = 1
analyze = 1
makefig = 1

diffuse = 0
antibunch = 1
pulsed = 1

#General saving folder parameters
filepath = "/mnt/c/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/DotTransferSim/"
filedir = "Sens-Tau-LP-100ligands/"
filenames = [""]

#file length - goes to the max of these two
numlines = 10**8
maxlines = 5*10**9
endtime = 3*10**14 #ps

#Sample parameters
temp = 298 #K
k_demissioni = 1400000 #emission lifetime in ps
k_sem = 10000 #short lived singlets
k_tem = 100000000 #long lived triplets
k_fiss = 1 #made up from marks paper on pentacene fission
k_trans = 200000 #based on Rao paper 2018
emwavelength = 815
r = 10 #nm - hydrodynamic radius of particles

eta = 8.9 * 10**(-13) # kg/nm s - dynamic viscosity of solvent (water)
n = 1.3 # index of refraction - sample in water

concentration = 2*10**(-10) 
nligands = 100

dabsXsec = 7.180616773321853*10**(-9)# per emitter numabs'd = phperpulse*absXsec*numEms - this is reasonable based on absXsec for CdSe is 550000*r^3/cm (from Bawendi paper Ruvim sent me)
labsXsec = dabsXsec*10

#Laser parameters
reprate = 0.05 #MHz

wavelength = 532 #nm
beamdiam = 5000000 #5 mm in nm
laserpwri = 0.5 # mW (0.52 mWinto back of objective 23 #mW)
#laserpwr = [0.05,0.1,0.2,0.3,0.5,0.75,1,10]#[0.001,0.01,0.025,

pulselength = 80 #ps - not used


#Objective parameters
foclen = 310000 #310 microns in nm (working distance + coverslip thickness)
NA = 1.4

#Detector parameters
darkcounts = 100 #s^-1
sensitivity = 0.01 
deadtime = 70#70000 #70 ns in ps
afterpulse = 0.0001 #percent of time a photon is emitted a deadtime after one is detected

#Correlation parameters
order = 2 #g2
mode = "t2"
gnpwr = 20
numbins = 4096
pulsebins = (10**2)-1#should always be an odd number
channels = 2

#Miscellaneous simulation parameters
picyzoom = [-1,-1]
timestep = 1000000 #average number of photons per "round" of calculations
for i in range(3):
    sensitivity = 0.01/8
    laserpwr = laserpwri*10**(i-1)
    for j in range(3):
        sensitivity = sensitivity*8
        for k in range(3):
            k_demission = k_demissioni*10**(k-1)
            
            filename = "sens-" + str(int(100*sensitivity)) + "-lp-"+str(int(1000*laserpwr)) + "uW-tauDotEm-"+str(int(k_demission/1000)) +"ns"

            sim.simulate(filepath, filedir, filename, write, analyze, makefig, diffuse, antibunch,
                            pulsed, numlines, maxlines, endtime,temp, concentration, 
                            dabsXsec, labsXsec, k_demission, k_sem,
                            emwavelength, r,eta, n, k_tem, k_fiss, k_trans, 
                            reprate, wavelength, laserpwr, pulselength, foclen,
                            NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, order, 
                            mode, gnpwr, numbins, pulsebins, channels)
