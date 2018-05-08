import scipy
import numpy

import diffuse as d

def newcphoton(antibunch, k_emission):
    if antibunch == 1:
        add = scipy.random.exponential(k_emission)
    #elif antibunch == 2:
    #    add = scipy.random.exponential(k_emission+k_excitation)
    else:
        add = - (k_emission*numpy.log(1-numpy.random.rand())) #poisson
        
    return add

def quicksort(array):
    less = []
    pivot = []
    more = []
    if len(array) <= 1:
        return array
    else:
        pivotpt = array[int(len(array)/2)]
        for i in array:
            if i < pivotpt:
                less.append(i)
            elif i > pivotpt:
                more.append(i)
            else:
                pivot.append(i)
        less = quicksort(less)
        more = quicksort(more)
    return less + pivot + more
#recursion depth limit is around 1000

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

def excitepulsed(k_ex, probex, nextem, taurep):
    return int(nextem/taurep)*taurep + numpy.random.geometric(p=probex)*taurep      

def excitects(k_ex, probex, nextem, taurep):
    return nextem + scipy.random.exponential(k_ex)

def ff(val, biggerval, otherval, lifetime, factor):
    while val < biggerval and val < otherval:
        if biggerval - val < factor*lifetime:#fast forward close to the right time
            val = int((biggerval - val)/lifetime)*lifetime + val
        val = val + scipy.random.exponential(lifetime)
    return val

def nextphotonss(lastphoton, sensitivity,
               k_demission, k_fiss, k_trans, k_sem, k_tem, k_dexcitation, k_lexcitation, 
               diffouttime, numEms, nextOut, nextIn, endround,
               antibunch, pulsed, taurep,
               dabsXsec, labsXsec, photonsperpulse, AvgEms, diffsIn, 
               diffsOut, ndiffsOut):
    
    nextphoton = []
    
    newphoton = newcphoton
    if pulsed == 1:
        excite = excitepulsed
    else:
        excite = excitects
    probdex = 1- (1-dabsXsec)**photonsperpulse #1-probability of not being excited each pulse
    problex = 1- (1-labsXsec)**photonsperpulse #1-probability of not being excited each pulse
     
     
    ''' 
    while nextOut < endround: #calculate emissions before it diffuses out
        #print("1")
        diffOut = numpy.random.randint(numEms)#index identity of emitter which diffused out
        nextem = lastphoton[diffOut] #lastphoton is an array of length numEms holding the time of last photon emitted for each emitter
        

        while nextem < nextOut:
            #First, excite:
            nextdex = excite(k_dexcitation, probdex, nextem, taurep)
            nextlex = excite(k_lexcitation, problex, nextem, taurep)
            
            #regular emission
            #nextem = nextem + newphoton(antibunch, k_demission)

            if nextdex < nextlex: # dot must emit before 
                nextem = nextdex + scipy.random.exponential(k_demission)
                if numpy.random.rand() < sensitivity and nextem < min(nextOut): # if we didn't miss it due to sensitivity
                    nextphoton = insert(nextphoton, nextem)

            fisstime = scipy.random.exponential(k_fiss)
            transtime = scipy.random.exponential(k_trans)
            semtime = scipy.random.exponential(k_sem)
            temtime = scipy.random.exponential(k_tem)

            if semtime < fisstime:
                nextem = nextlex + semtime
                if numpy.random.rand() < sensitivity and nextem < min(nextOut): # if we didn't miss it due to sensitivity
                    nextphoton = insert(nextphoton, nextem)
            
            else:
                if transtime < temtime: # first one didn't emit before transfer - we assume (shortpass) filter triplet emission so we only see dot photons - regardless, triplet emission should basically always be slower than transfer
                    nextem = nextlex + transtime + scipy.random.exponential(k_demission)

                    if numpy.random.rand() < sensitivity and nextem < min(nextOut): # if we didn't miss it due to sensitivity
                        nextphoton = insert(nextphoton, nextem)
                    
                transtime = scipy.random.exponential(k_trans)
                temtime = scipy.random.exponential(k_tem)
                if transtime < temtime: # second one didn't emit before transfer
                    nextem = nextlex + transtime + scipy.random.exponential(k_demission)

                    if numpy.random.rand() < sensitivity and nextem < min(nextOut): # if we didn't miss it due to sensitivity
                        nextphoton = insert(nextphoton, nextem)
                

        numEms = numEms - 1
        del lastphoton[diffOut]
        diffsOut = diffsOut + 1
        nextOut = nextOut + d.diffuse(diffouttime, numEms)
        if numEms == 0:
            break #we'll want to start a new round
        
    while nextIn < endround:
        lastphoton.append(nextIn)
        diffsIn = diffsIn + 1
        numEms = numEms + 1
        nextIn = nextIn + d.diffuse(diffouttime, AvgEms)
        if numEms == 1:
            endround = nextIn
            break #skip to next round if we just diffused out and back in'''


    #Note that as it stands, this only works for: 
    # 1) STATIONARY (non-diffusing) - could be fixed by adding in the possibility of the 
    # dot diffusing out before emission. Just need to be careful about the timings everywhere 
    # and so there might be too many "if"s slowing down the code. Also need to un-comment
    # out the diffusion-related bits of code.
    # 
    # 2) PULSED - could be fixed by allowing the dot to be re-excited after emitting 
    # before another triplet transfer. We would also have to consider new ligand excitations.
    # 
    # 3) SINGLE DIMER - could be fixed by building in the possibility of multiple
    # ligands trying to transfer to the same dot. This might be horribly slow.     
    for i in range(numEms):
        #print(2)
        nextem = lastphoton[i]

        #print(3)
        while nextem < endround: #continue adding photons until end of round, on average 10 emissions
            
            fisstime = 0
            semtime = -1
            temtime = -1
            transtime = 0
            dotem = 1
            transtime = [1,1]
            temtime = [0,0]
            trans2time = 1
            tem2time = 0
            #assumes triplet lifetime is much longer than dot lifetime or transfer time 
            while (nextdex < min(transtime) and dotem > max(temtime)) or (transtime[0] > temtime[0] and transtime[1] > temtime[1]):
                #first successfully get to a fissed state:
                fisstime = semtime + 1
                while semtime < fisstime:
                    #if we're in a later round and the dot emitted before ligand fission, we'll see that
                    if dotem < semtime and numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                        nextphoton = insert(nextphoton, dotem)
                        nextem = dotem

                        
                    #First, excite:
                    else: #re-excite the dot and find it's emission time
                        nextdex = excite(k_dexcitation, probdex, nextem, taurep)
                        dotem = nextdex + scipy.random.exponential(k_demission)

                    #re-excite the ligands
                    nextlex = excite(k_lexcitation, problex, nextem, taurep)
                    #and calculate the singlet emission time and the fission time
                    semtime = nextlex + scipy.random.exponential(k_sem)
                    fisstime = nextlex + scipy.random.exponential(k_fiss)

                #now we've successfully done fission (the dot may have emitted in the interim and we know when it emits next)

                transtime= [fisstime + scipy.random.exponential(k_trans),fisstime + scipy.random.exponential(k_trans)]
                temtime = [fisstime + scipy.random.exponential(k_tem),fisstime + scipy.random.exponential(k_tem)]

            #now atleast one triplet successfully transfers to the dot to emit

            index = transtime.index(min(transtime))
            otherone = (index + 1)%2

            if nextdex < transtime[index]:
                #dot emits, no transfers until it does
                if numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                    nextphoton = insert(nextphoton, dotem)
                    nextem = dotem
                if transtime[index] > nextem:
                    nextem = transtime[index] + scipy.random.exponential(k_demission)#excitation transfers to dot and emits
                    if numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                        nextphoton = insert(nextphoton, nextem)
                    
                    #fast forward close to the right time
                    transtime[otherone] = ff(transtime[otherone], nextem, temtime[otherone], k_trans, 5)
                    
                    
                    if transtime[otherone] < temtime[otherone]:
                        nextem = transtime[otherone] + scipy.random.exponential(k_demission) #excitation transfers to dot and emits
                        if numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                            nextphoton = insert(nextphoton, nextem)
                else: # the first transfer "happens" (tries to) while the dot is still excited
                    #fast forward close to the right time
                    transtime[otherone] = ff(transtime[otherone], nextem, temtime[otherone], k_trans, 5)
                    transtime[index] = ff(transtime[index], nextem, temtime[index], k_trans, 5)
                    #now both transfer times are bigger than nextem OR less than their temtime

                    #if we still have at least one of them that transfers
                    if not (transtime[index] > temtime[index] and transtime[otherone] > temtime[otherone]):
                        index = transtime.index(min(transtime))
                        otherone = (index + 1)%2
                        #the first one emits
                        nextem = transtime[index] + scipy.random.exponential(k_demission)#excitation transfers to dot and emits
                        if numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                            nextphoton = insert(nextphoton, nextem)
                        #the second one has a chance to emit
                        transtime[otherone] = ff(transtime[otherone], nextem, temtime[otherone], k_trans, 5)

                        if transtime[otherone] < temtime[otherone]:
                            nextem = transtime[otherone] + scipy.random.exponential(k_demission) #excitation transfers to dot and emits
                            if numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                                nextphoton = insert(nextphoton, nextem)
            
            #OTHERWISE, the triplet transfer occurred before the dot could get excited
            else:
                nextem = transtime[index] + scipy.random.exponential(k_demission)#excitation transfers to dot and emits
                if numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                    nextphoton = insert(nextphoton, nextem)
                #the second one has a chance to emit
                transtime[otherone] = ff(transtime[otherone], nextem, temtime[otherone], k_trans, 5)

                if transtime[otherone] < temtime[otherone]:
                    nextem = transtime[otherone] + scipy.random.exponential(k_demission) #excitation transfers to dot and emits
                    if numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                        nextphoton = insert(nextphoton, nextem)

        lastphoton[i] = nextem
        
        
'''
    if numEms == 0:
        numEms = 1
        while nextOut<nextIn:
            nextOut = nextOut + d.diffuse(diffouttime, 1)
            ndiffsOut = ndiffsOut + 1
        lastphoton = [nextIn]
        nextIn = nextIn + d.diffuse(diffouttime, AvgEms)
        diffsIn = diffsIn + 1'''

    return nextphoton, lastphoton, numEms, nextIn, nextOut, diffsIn, diffsOut, ndiffsOut, 

def darkcounts(avtime, length, lastdc, deadtime):
    dc = []
    while len(dc) < length:
        new = lastdc -((avtime)/2)*(numpy.log((1-numpy.random.rand())))#so that the interval is (0,1] instead of [0,1)
        #divide avtime by 2 because assuming 2 detectors
        if not new - lastdc < deadtime:
            dc.append(new)
        elif numpy.random.rand() < 0.5:#if the darkcounts appeared within the deadtime of each other but on different detectors
            dc.append(new)
        lastdc = new
    return dc
