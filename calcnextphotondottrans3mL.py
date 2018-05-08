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

def fastforward(transtime, nextem, temtime, k_trans, factor):
    for i in range(len(transtime)):
        for j in range(2):
            transtime[i][j] = ff(transtime[i][j], nextem, temtime[i][j], k_trans, factor)
    return transtime


def alltransgttemtime(transtime, temtime):
    for i in range(len(transtime)):
        if transtime[i][0] < temtime[i][0] or transtime[i][1] < temtime[i][1]:
            return 0
    return 1

def deltransgttem(transtime,temtime):
    d = 0
    for i in range(len(transtime)):
        if transtime[i-d][0] > temtime[i-d][0]:
            if transtime[i-d][1] > temtime[i-d][1]:
                del transtime[i-d]
                del temtime[i-d]
                d = d + 1
            else:
                transtime[i-d][0] = float("inf")
        elif transtime[i-d][1] > temtime[i-d][1]:
            transtime[i-d][1] = float("inf")
    return transtime, temtime


def allsemltallfiss(fisstime,semtime):
    for i in range(len(fisstime)):
        if fisstime[i] < semtime[i]:
            return 0               
    return 1

def findmin(array):
    minval = array[0][0]
    index = [0,0]
    for i in range(len(array)):
        if array[i][0] < minval:
            minval = array[i][0]
            index = [i,0]
        if array[i][1] < minval:
            minval = array[i][1]
            index = [i,1]
    return index, minval

def findmax(array):
    minval = array[0][0]
    index = [0,0]
    for i in range(len(array)):
        if array[i][0] > minval:
            minval = array[i][0]
            index = [i,0]
        if array[i][1] > minval:
            minval = array[i][1]
            index = [i,1]
    return index, minval

def nextphotonss(lastphoton, sensitivity, nligands,
               k_demission, k_fiss, k_trans, k_sem, k_tem, k_dexcitation, k_lexcitation, 
               diffouttime, numEms, nextOut, nextIn, endround,
               antibunch, pulsed, taurep,
               dabsXsec, labsXsec, photonsperpulse, AvgEms, diffsIn, 
               diffsOut, ndiffsOut, testdummy):
    
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
  
    for j in range(numEms):
        nextem = lastphoton[j]
        #print(3)
        while nextem < endround: #continue adding photons until end of round, on average 10 emissions
            transtime = []
            temtime = []
            fisstime = []
            semtime = []
            for i in range(nligands):
                nextfisstime = 1
                nextsemtime = 0
                while nextsemtime < nextfisstime:
                    nextlex = excite(k_lexcitation, problex, nextem + semtime[i], taurep)
                    nextfisstime = nextlex + scipy.random.exponential(k_fiss)
                    nextsemtime = nextlex + scipy.random.exponential(k_sem)
                fisstime.append(nextfisstime)
                semtime.append(nextsemtime)
                
                transtime.append([fisstime[i] + scipy.random.exponential(k_trans),fisstime[i] + scipy.random.exponential(k_trans)])
                temtime.append([fisstime[i] + scipy.random.exponential(k_tem),fisstime[i] + scipy.random.exponential(k_tem)])
                
            nextdex = excite(k_dexcitation, probdex, nextem, taurep)
            dotem = nextdex + scipy.random.exponential(k_demission)
            atgttt = alltransgttemtime(transtime, temtime)
            asltaf = allsemltallfiss(fisstime, semtime)
            index, minval = findmin(transtime)
            tindex, maxval = findmax(temtime)
            #assumes triplet lifetime is much longer than dot lifetime or transfer time 
            while (nextdex < minval and dotem > maxval) or atgttt == 1:
                #if testdummy > 0 :
                #print("first successfully get to a fissed state:" + str(testdummy))
                for i in range(nligands):
                    nextlex = excite(k_lexcitation, problex, nextem, taurep)
                    semtime[i] = nextlex + scipy.random.exponential(k_sem)
                    fisstime[i] = nextlex + scipy.random.exponential(k_fiss)
                #print(semtime)
                asltaf = allsemltallfiss(fisstime,semtime)

                while asltaf == 1:
                    #print("if we're in a later round and the dot emitted before ligand fission, we'll see that")
                    '''print("dotem < semtime = write:")
                    print(dotem)
                    print(semtime)
                    print("-----------------")'''
                    if dotem < min(semtime): 
                        nextem = dotem
                        if numpy.random.rand() < sensitivity:
                            nextphoton = insert(nextphoton, dotem)
                        nextdex = excite(k_dexcitation, probdex, nextem, taurep)
                        dotem = nextdex + scipy.random.exponential(k_demission)


                        

                    #re-excite the ligands
                    
                    #and calculate the singlet emission time and the fission time
                    for i in range(nligands):
                        nextlex = excite(k_lexcitation, problex, nextem, taurep)
                        semtime[i] = nextlex + scipy.random.exponential(k_sem)
                        fisstime[i] = nextlex + scipy.random.exponential(k_fiss)
                    #print(semtime)
                    asltaf = allsemltallfiss(fisstime,semtime)
                    

                #now we've successfully done fission (the dot may have emitted in the interim and we know when it emits next)

                for i in range(nligands):
                    if not semtime[i] < fisstime[i]:
                        transtime[i] = [fisstime[i] + scipy.random.exponential(k_trans),fisstime[i] + scipy.random.exponential(k_trans)]
                        temtime[i] = [fisstime[i] + scipy.random.exponential(k_tem),fisstime[i] + scipy.random.exponential(k_tem)]
                    else:
                        transtime[i] = float("inf")
                
                atgttt = alltransgttemtime(transtime, temtime)
                index, minval = findmin(transtime)
                tindex, maxval = findmax(temtime)
                #print(minval)
            #print(transtime)
            transtime, temtime = deltransgttem(transtime,temtime)
            '''if not transtime == [] #should never be empty because to get out of previous while loop, at least one temtime had to be longer than the transtime
                index, minval = findmin(transtime)
            else:
                minval = nextdex + 1
            if nextdex < minval:
                #if testdummy >0 :
                #    print("in the if 209" + str(testdummy))
                #dot emits, no transfers until it does
                if numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                    nextphoton = insert(nextphoton, dotem)
                    nextem = dotem'''

            #now atleast one triplet successfully transfers to the dot to emit, we hope
            while not transtime == []:
                #print("transtime isn't empty")
                print(transtime)
                index, minval = findmin(transtime)
                if nextdex < minval and nextdex > nextem:
                    if numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                        nextphoton = insert(nextphoton, dotem)
                        nextem = dotem
                    transtime = fastforward(transtime, nextem, temtime, k_trans, 5)
                    transtime, temtime = deltransgttem(transtime,temtime)
                    index, minval = findmin(transtime)

                if minval > nextem:
                    nextem = minval + scipy.random.exponential(k_demission)#excitation transfers to dot and emits
                    if numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                        nextphoton = insert(nextphoton, nextem)
                
                    transtime[index[0]][index[1]] = float("inf")
                
                transtime = fastforward(transtime, nextem, temtime, k_trans, 5)
                transtime, temtime = deltransgttem(transtime,temtime)
            print(nextem)#nextem is shrinking?



        lastphoton[j] = nextem
        
        
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
