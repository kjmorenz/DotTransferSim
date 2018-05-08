
import diffuse as d
import calcnextphotondottrans6 as calc1
import calcnextphotondottrans6mL as mLcalc
import time as rt
import math
import numpy

def insert(array, val):
    if array == []:
        return [val]
    if array[-1] < val:
        array.append(val)
        return array
    if array[0] > val:
        array.insert(0, val)
        return array
    index = int(len(array)/2)
    while array[index] < val:
        index = index + 1
    while array[index - 1] > val:
        index = index - 1
    array.insert(index, val)
    return array

def addafterpulse(photons, darks, deadtime, afterpulse):
    for i in range(len(photons)):
        apphoton = photons[i]+ deadtime
        index = 0
        if numpy.random.rand() < afterpulse:#includes probability that went to other detector
            if apphoton < photons[-1]:  #ignore afterpulsing after end of round
                while photons[i+index] < apphoton:
                    index = index + 1
                index = index - 1
                if apphoton - photons[index] > deadtime: #so long as the second detector wasn't already turned off by another photon arrival
                    photons.insert(index-1, apphoton)
            else:#add it to dark counts instead so it'll get added in next round
                darks = insert(darks, apphoton)
    return photons, darks



def write(filepath, filedir, fullfilename, antibunch, diffuse, pulsed, numlines, maxlines, endtime,
            temp, concentration, dabsXsec, labsXsec,k_demission,
            k_fiss, k_trans, k_sem, k_tem, emwavelength, r,
            eta, n, reprate,wavelength, laserpwr, pulselength, foclen,
            NA, darkcounts, sensitivity, nligands, deadtime, afterpulse, timestep, channels, seq, mL1):

    if nligands != 1 or mL1 == 1:
        calc = mLcalc
    else:
        calc = calc1


    h = 6.62607004*10**(-34) # m^2 kg / s (Planck's const)
    c = 299792458 # m/s (Speed of light)
    N_A = 6.023*(10**23) # Avogadro
    K_B = 1.38064852*10**(-5) # Boltzmann in kg nm^2 /K s^2
    suffix = ".txt"
    #beamwaist = 2*wavelength*foclen/(math.pi*beamdiam) #nm
    #this gives too small a value (~20 nm) - gets us very fast diffusion times, too fast to be real
    ##beamwaist = 4lambda f / 2 pi D where f is the focal length and D is the diameter of light entering
    ## or lambda/2NA

    beamwaist = wavelength/(2*NA) #~200 nm

    focalVol = (math.pi**(1.5))*(2.33*n/NA)*beamwaist**3 #nm^3
    AvgEms = (concentration * N_A * focalVol / (10**24)) #number not mol
    taurep = (10**12)/(reprate*10**6) #ps

    energyperphoton = h*c/(wavelength/(10**9)) #J


    energyperpulse = laserpwr/(reprate * 1000000000) #J
    phperpulse = energyperpulse/energyperphoton
    phpers = laserpwr/(energyperphoton*1000)


    #make rounds average timestep emissions per round:
    if pulsed == 1:
        k_dexcitation = 10**6/(phperpulse*dabsXsec*reprate)#*AvgEms) #ps/excitation
        k_lexcitation = 10**6/(phperpulse*labsXsec*reprate)#*AvgEms) #ps/excitation
        experemperpulse = dabsXsec*phperpulse
        exsperpulse = dabsXsec*phperpulse*AvgEms + labsXsec*phperpulse*nligands*AvgEms
        AvgExEvs = exsperpulse*reprate*10**6
        if experemperpulse > 1:
            experemperpulse = 1
        experps = experemperpulse*taurep*AvgEms # comment *avgems out to make per em
        timestep = timestep/experps #ps per round
        if timestep < taurep:
            timestep = taurep

    #CW version:
    else:
        k_dexcitation = 10**12/(phpers*dabsXsec)#*AvgEms) #ps/excitation
        k_lexcitation = 10**12/(phpers*labsXsec)#*AvgEms) #ps/excitation
        AvgEmEvs = AvgEms/(k_demission + k_dexcitation) # average sample emissions per second (in pHz)
        AvgExEvs = AvgEms*phpers*(dabsXsec+labsXsec*nligands)
        timestep = timestep/AvgEmEvs # in ps
    print("kdex = " + str(k_dexcitation))
    print("klex = " + str(k_lexcitation))
    endround = timestep
    print(timestep)


    ##focal vol = pi^1.5 kw^3, k is a geometric factor usually taken as optical resolution in the z direction
    ##divided by that in the xy direction, w is the e^-2 beam radius in the xy plane
    ##k >= 2.33 n/NA, r_z = 2nlambda/NA^2, rxy = 0.61 lambda/NA
    ##http://www.fcsxpert.com/classroom/theory/what-is-confocal-volume.html

    diffouttime = (3*((beamwaist)**2)*(math.pi**2)*eta*r/(8*K_B*temp))*10**12
    #= 3 (b pi)^2 eta r / 8 k_B T
    ##sigma = numpy.sqrt(2*beamwaist/(6*math.pi*eta*r))

    #gets a similar order of magnitude if x = focalVol^(2/3)

    ##Phys Chem by Atkins pg 771 - prob of being a distance x from origin after time t is
    ##P = sqrt(2tau/pi t)*e^(-x^2 tau/(2t lambda^2) where tau is the time it takes to go a distance lambda
    ##Diffusion coeff D = lambda^2/2tau = uRT/zF = uK_BT where u is the mobility or the ratio of the
    ##terminal drift velocity to an applied force
    ##lambda = zuF, also u = 1/ 6 pi eta r Stokes Einstein
    ##
    ##Better: Wikipedia, (also in Atkins) Einstein relation:
    ## D = k_B T/(6 pi eta r), eta is the viscosity, r is the radius of the spherical particle
    ##
    ##or pg 770, mean distanced travelled by particles of radius r in solvent of viscosity eta is
    ##                       x = sqrt(2kTt/3 pi^2 eta r)

    ##https://www2.gwu.edu/~phy21bio/Reading/randomwalkBerg.pdf gives sigma_x is sqrt(2Dt), P(x) = (sqrt(4 pi Dt)^-1)e^((-x^2)/4Dt))

    # combination of Atkins and Wikipedia gives:
    # mean free path lambda = (root(2)*conc*collision cross sectional area)
    # --> (using equations about D above) Mean free time = 3 eta / 2 conc ^2 pi r^3 Kb T where conc is only acc dependent
    # and r is the sum of the two diameters


    deltat = 0
    line = 0
    sigcts = 0
    nextOut = 0
    nextIn = 0
    diffsIn = 0
    diffsOut = 0
    ndiffsOut = 0
    lastwritten = -deadtime - 1 #both detectors on
    lastlastwritten = lastwritten
    dclen = 20
    channel = numpy.random.randint(channels)

    testdummy =0



    numEms = round(AvgEms) #start with the average number of emitters in the focal volume
    if numEms <= 0: #always start with at least one emitter in the focal volume
        numEms = 1
    print("NumEms = " + str(numEms))
    print("AvgEms = " + str(AvgEms))



    lastphoton =[]
    for i in range(numEms):
        lastphoton.append(0)

    nextdex = numpy.random.geometric(p=1- (1-dabsXsec)**phperpulse)*taurep

    lastdiff = 0 # for calculating actual avg ems
    tnumEms = 0

    if diffuse == 1:
        nextOut = d.diffuse(diffouttime, numEms)#the number which have a chance to diffuse out are the number in
        nextIn = d.diffuse(diffouttime, AvgEms)#the number which have a chance to diffuse in are the average number
        nextdiff = min(nextIn, nextOut)

    elif diffuse == 0:
        nextdiff = float('inf')
        nextOut = float('inf')
        nextIn = float('inf')

    lastsyncpulse = 0
    if pulsed != 1:
        lastsyncpulse = float('inf')

    avtime = (1/darkcounts)*(10**12) #average time between darkcount photons in ps given darkcounts in s^-1
    darks = calc.darkcounts(avtime, dclen, 0, deadtime)

    dcpointer = 0

    file = open(filepath +"RawData/"+ filedir+ fullfilename +"/" + fullfilename+ suffix, 'w')
    difffile = open(filepath +"RawData/"+ filedir+ fullfilename +"/" + "diffsof" + fullfilename+ suffix, 'w')
    starttime = rt.time()
    wrotediff = [0,0]
    while (line < numlines or lastwritten < endtime) and line < maxlines:
        if line > maxlines:
            break
        photons = []



        if dcpointer == dclen:
            darks = calc.darkcounts(avtime, dclen, darks[dcpointer-1], deadtime)
            dcpointer = 0

        if numEms < 1:
            print("NumEms = " + str(numEms))
            print("AvgEms = " + str(AvgEms))
            print("Rounded = " + str(round(AvgEms)))
            please = crash


        if len(lastphoton) != numEms :
            print("numEms:" + str(numEms))
            print("lastphoton:" +str(lastphoton))
            print("write line 206")

        if wrotediff != [nextIn, nextOut] and nextIn < nextOut:
            difffile.write(str(nextIn) + ","+ str(numEms + 1)+"\n")
            wrotediff = [nextIn, nextOut]
        elif wrotediff != [nextIn, nextOut] and nextOut < nextIn and numEms >=0:
            difffile.write(str(nextOut) + ","+ str(numEms - 1)+"\n")
            if numEms -1 == 0:
                difffile.write(str(nextIn) + ",1\n")
            wrotediff = [nextIn, nextOut]
        #print("New Round")
        #print("sensitivity is " + str(sensitivity))
        dataphotons, lastphoton, numEms, nextIn, nextOut, diffsIn, diffsOut, ndiffsOut, nextdex = calc.nextphotonss(lastphoton, sensitivity, nligands,
                                                                        k_demission, k_fiss, k_trans, k_sem, k_tem, k_dexcitation, k_lexcitation,
                                                                        diffouttime, numEms, nextOut, nextIn, endround,
                                                                        antibunch, pulsed, taurep,
                                                                        dabsXsec, labsXsec, phperpulse, AvgEms, diffsIn,
                                                                        diffsOut, ndiffsOut, testdummy, nextdex,seq)
        #if not dataphotons == []:
        #print(dataphotons)
        endround = endround + timestep
        if len(lastphoton) == 1 and endround < lastphoton[0]:
            endround = lastphoton[0]
    ##            print("new dataphotons")
    ##            print(dataphotons)


        if not dataphotons == []:
            sigcts += len(dataphotons)

            datapointer = 0

            while datapointer < len(dataphotons):
                if dcpointer == dclen:
                    darks = calc.darkcounts(avtime, dclen, darks[dcpointer-1], deadtime)
                    dcpointer = 0

                if darks[dcpointer] < dataphotons[datapointer]:
    ##                        print(0)
                    if darks[dcpointer] - lastwritten > deadtime or (darks[dcpointer] - lastlastwritten > deadtime and numpy.random.rand()>0.5):
                        photons.append(darks[dcpointer])
                        lastlastwritten = lastwritten
                        lastwritten = darks[dcpointer]
    ##                            print("wrote a dc - 3")

                    dcpointer = dcpointer + 1


                else:
    ##                            if datapointer >= len(dataphotons) or photons == []:
    ##                                print(darks)
    ##                                print(dcpointer)
    ##                                print(dataphotons)
    ##                                print(datapointer)
    ##                                print(photons)
    ##                                print(dataphotons[datapointer] + photons[0])
                    if dataphotons[datapointer] - lastwritten > deadtime or (dataphotons[datapointer] - lastlastwritten > deadtime and numpy.random.rand()>0.5):
                        photons.append(dataphotons[datapointer])
                        lastlastwritten = lastwritten
                        lastwritten = dataphotons[datapointer]
    ##                            print("wrote data - 1")

                    datapointer = datapointer + 1
    ##                        print(datapointer)




        #the first time through, we calculate how often to print progress
        if deltat == 0:
            deltat = (rt.time() - starttime)*1000
            print("Delta t for one round: " + str(deltat))
            print(timestep)
            #print(photons)

        if not deltat == 0:
            testdummy = testdummy  + 1
        if not photons == []:
            if not afterpulse == 0:
                photons, darks = addafterpulse(photons, darks, deadtime, afterpulse)
            #writing into file
            if photons[0] < 0:
                print(line + " is neg - didn't write")
            elif photons[0] == float('inf'):
                print(line + " is inf - didn't write")
            else:
                while photons[0] > lastsyncpulse:
                    file.write(str(channels) + "," + str(lastsyncpulse)+"\n")
                    lastsyncpulse = lastsyncpulse + int(taurep)
                    line = line + 1

                if photons[0] - lastwritten <= deadtime: #leq since afterpulse photons have to go to other det too
                    channel = (channel+1)%channels #went to the other detector
                else:
                    channel = numpy.random.randint(channels)
                file.write(str(channel) + "," + str(int(photons[0]))+"\n")
    ##            print("wrote a line")
                line = line + 1
                if line%numlines==0:
                    print("time = " + str(lastwritten/10**12) + " s")
                    if lastwritten>endtime:
                        break #should never happen
                if line%(numlines/10) == 0:
                    print(fullfilename + ": " + str(line*100/numlines) + "% - " + str(photons[0]/10**9) + " ms")
                elif deltat > 1 and line%(numlines/100)==0:#if it's slow print more often
                    print(str(line*100/numlines) + "% - " + str(photons[0]/10**9) + " ms")


            for j in range(1,len(photons)):
                if photons[j] < 0:
                    print(line + " is neg - didn't write")
                elif photons[j] == float('inf'):
                    print(line + " is inf - didn't write")
                else:
                    while photons[j] > lastsyncpulse:
                        file.write(str(channels) + "," + str(lastsyncpulse)+"\n")
                        lastsyncpulse = lastsyncpulse + int(taurep)
                        line = line + 1
                        if line%numlines==0:
                            print("time = " + str(lastwritten/10**12) + " s")
                            if lastwritten>endtime:
                                break #should never happen
                        if line%(numlines/10) == 0:
                            print(fullfilename + ": " + str(line*100/numlines) + "% - " + str(photons[j]/10**9) + " ms")
                        elif deltat > 1 and line%(numlines/100)==0:#if it's slow print more often
                            print(str(line*100/numlines) + "% - " + str(photons[j]/10**9) + " ms")

                    if photons[j] - photons[j-1] <= deadtime: #leq since afterpulse photons have to go to other det too
                        channel = (channel+1)%channels #went to the other detector
                    else:
                        channel = numpy.random.randint(channels)
                    file.write(str(channel) + "," + str(int(photons[j]))+"\n")
        ##            print("wrote a line")
                    line = line + 1
                    if line%numlines==0:
                        print("time = " + str(lastwritten/10**12) + " s")
                        if lastwritten>endtime:
                            break #should never happen
                    if line%(numlines/10) == 0:
                        print(fullfilename + ": " + str(line*100/numlines) + "% - " + str(photons[j]/10**9) + " ms")
                    elif deltat > 1 and line%(numlines/100)==0:#if it's slow print more often
                        print(str(line*100/numlines) + "% - " + str(photons[j]/10**9) + " ms")

            #print("Lines written: " + str(line))






    print("Last photon time = " + str(int(max(lastphoton))/(10**12)) + " s")
    print("Signal counts = " + str(sigcts))
    print("Diffusion-in events = " + str(diffsIn))
    print("Diffusion-out events = " + str(diffsOut))
    print("Final emitters in = " + str(numEms))
    file.close()
    difffile.close()
    file = open(filepath +"Figures/"+ filedir + fullfilename +"/" + "params"+ suffix, 'w')
    file.write(fullfilename +" - date and time: " + str(rt.time()) + "\n")
    file.write("Last photon time = " + str(int(max(lastphoton))/(10**12)) + " s"+"\n")
    file.write("Signal counts = " + str(sigcts)+"\n")
    file.write("Excitations per pulse " + str(exsperpulse) +"\n")
    file.write("Excitations per sec " + str(AvgExEvs) +"\n")
    file.write("Diffusion-in events = " + str(diffsIn)+"\n")
    file.write("Diffusion-out events = " + str(diffsOut)+ ", " +str(ndiffsOut) + "\n")
    file.write("Final emitters in = " + str(numEms) + "\n")
    file.write("Intended Avg Ems = " + str(AvgEms) + "\n")
    file.write("Number of records = " + str(numlines)+"\n")
    file.write("k_demission = " + str(int(k_demission))+"\n")
    file.write("k_dexcitation = " +str(int(k_dexcitation))+"\n")
    file.write("dabsXsec = " + str(dabsXsec) + "\n")
    file.write("Ligands per dot = " +str(nligands) +"\n")
    file.write("k_lemission = " + str(int(k_sem))+"\n")
    file.write("k_lexcitation = " +str(int(k_lexcitation))+"\n")
    file.write("labsXsec = " + str(labsXsec) + "\n")
    file.write("darkcounts = " + str(darkcounts) + " Hz \n")
    file.write("sensitivity = " + str(sensitivity) + "\n")
    file.write("deadtime = " + str(deadtime/1000) + " ns \n")
    file.write("afterpulsing = " + str(afterpulse) +"\n")
    file.write("radius = " + str(r)+" M\n")
    file.write("concentration = " + str(concentration)+" M\n")
    file.write("Solvent viscosity = " + str(eta)+" \n")
    file.write("Solvent refractive index = " + str(n)+" \n")
    file.write("Dark fission lifetime = " + str(k_fiss)+" \n")
    file.write("Triplet emission lifetime = " + str(k_tem)+" \n")
    file.write("Transfer lifetime = " + str(k_trans)+" \n")
    file.write("Focal length = " + str(foclen)+" \n")
    file.write("NA = " + str(NA)+" \n")
    file.write("Laser power = " + str(laserpwr) + "\n")
    file.write("Time per round: " + str(timestep) + " ps \n")
##    file.write("Time limits for photon_gn = " + str(2**gnpwr) + "ps \n")
##    file.write("Pulse bins = " +str(pulsebins) + "\n")
    file.write("Average time to diffuse out of " + str(focalVol) + " nm^3 focal volume is " + str(diffouttime/(10**9)) + " \n")
    file.write("Antibunch? " + str(antibunch)+"\n")
    file.write("Diffusing as t = -diffouttime*numpy.log(1-numpy.random.rand()) #poissonian with average time diffouttime \n")
    file.write("\n")
    file.write("Main file + definitive input parameters:")
    # mf = open("C:/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/DotTransferSim/main.py", 'r')
    # for line in mf.readlines():
    #     file.write(line)
    # mf.close()
    file.close()
