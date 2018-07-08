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

def simulquicksort(array, b = [], c = [], d = [], e = []):
    less = []
    pivot = []
    more = []

    bless = []
    bpivot = []
    bmore = []

    cless = []
    cpivot = []
    cmore = []

    dless = []
    dpivot = []
    dmore = []

    eless = []
    epivot = []
    emore = []

    fless = []
    fpivot = []
    fmore = []
    if len(array) <= 1:
        return array, b, c, d, e
    else:
        pivotpt = array[int(len(array)/2)]
        for i in range(len(array)):
            val = array[i]
            if val < pivotpt:
                less.append(val)
                bless.append(b[i])
                cless.append(c[i])
                dless.append(d[i])
                eless.append(e[i])
            elif val > pivotpt:
                more.append(val)
                bmore.append(b[i])
                cmore.append(c[i])
                dmore.append(d[i])
                emore.append(e[i])
            else:
                pivot.append(val)
                bpivot.append(b[i])
                cpivot.append(c[i])
                dpivot.append(d[i])
                epivot.append(e[i])
        less, bless, cless, dless, eless = simulquicksort(less, bless, cless, dless, eless)
        more, bmore, cmore, dmore, emore = simulquicksort(more, bmore, cmore, dmore, emore)
    return less + pivot + more, bless + bpivot + bmore, cless + cpivot + cmore, dless + dpivot + dmore, eless + epivot + emore
    
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
       

def nextphotonss(lastphoton, sensitivity,
               k_emission, k_excitation, diffouttime,
               numEms, nextOut, nextIn, endround,
               antibunch, pulsed, taurep,bfrate,
               absXsec, photonsperpulse, AvgEms, diffsIn, 
               diffsOut, ndiffsOut, countbright, dftime = float("inf"), tripem = 0, 
               crashtime = float("inf"), crashtransferprob = 1, lastcrash = []):
    
    nextphoton = []
    
    newphoton = newcphoton
    probex = 1- (1-absXsec)**photonsperpulse #1-probability of not being excited each pulse
     
     
    while nextOut < endround: #calculate emissions before it diffuses out
        #print("1")
        diffOut = numpy.random.randint(numEms)#index identity of emitter which diffused out
        nextem = lastphoton[diffOut] #lastphoton is an array of length numEms holding the time of last photon emitted for each emitter
        

        while nextem < nextOut:
            #First, excite:
            if pulsed == 1:
                nextdex = int(nextem/taurep)*taurep + numpy.random.geometric(p=probdex)*taurep 
                nextlexs = []
                for ligands in numligands:
                    nextlexs = insert(nextlexs, int(nextem/taurep)*taurep + numpy.random.geometric(p=problex)*taurep)
                #geometric distribution gives us number of pulses it took to get excited
            else:
                nextdex = nextem + scipy.random.exponential(k_dexcitation)
                nextlexs = []
                for ligands in numligands:
                    nextlexs = insert(nextlexs, nextem + scipy.random.exponential(k_lexcitation))
            
            #regular emission
            #nextem = nextem + newphoton(antibunch, k_demission)

            if nextdex < min(nextlexs): # dot must emit before 
                nextem = nextdex + scipy.random.exponential(k_demission)
                if numpy.random.rand() < sensitivity and nextem < min(nextOut): # if we didn't miss it due to sensitivity
                    nextphoton = insert(nextphoton, nextem)
            
            fisstime = []
            transtime = []
            temtime = []

            index = 0
            for i in range(len(nextlexs)):
                while index < len(nextlexs):
                    nextlex = nexlexs[index]
                    thisfisstime = nextlex + scipy.random.exponential(k_fiss)
                    thissemtime = nextlex + scipy.random.exponential(k_sem)
                    
                    if thissemtime < thisfisstime:
                        ph = nextlex + semtime
                        if numpy.random.rand() < sensitivity and ph < min(nextOut): # if we didn't miss it due to sensitivity
                            nextphoton = insert(nextphoton, ph)
                        del nextlexs[index] 
                        index = index - 1
                    else:
                        transtime.append(nextlex + thisfisstime + scipy.random.exponential(k_trans))
                        fisstime.append(thisfisstime)
                        temtime.append(nextlex + thisfisstime + scipy.random.exponential(k_tem))
                    index = index + 1
            
            transtime, fisstime, temtime, nextlexs = simulquicksort(transtime, fisstime, temtime, nextlexs) #sort them now by the first to transfer to the dot
            
            for i in range(len(nextlexs)):   
                if transtime[i] < temtime[i]: # first one didn't emit before transfer - we assume (shortpass) filter triplet emission so we only see dot photons - regardless, triplet emission should basically always be slower than transfer
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
        lastcrash.append(nextIn)
        diffsIn = diffsIn + 1
        numEms = numEms + 1
        nextIn = nextIn + d.diffuse(diffouttime, AvgEms)
        if numEms == 1:
            endround = nextIn
            break #skip to next round if we just diffused out and back in
        
    for i in range(numEms):
        #print(2)
        nextem = lastphoton[i]
        bfem = 0
        
        if crashtime != float("inf"):
            nextcrash = lastcrash[i] 
        '''else:
            print("Error: crashtime is infinite in dda")
            crash = please '''

        #print(3)
        while nextem < endround: #continue adding photons until end of round, on average 10 emissions
            nextem = max(bfem, nextem)
            bfem = 0
            #First, excite:
            if pulsed == 1:
                nextem = int(nextem/taurep)*taurep + numpy.random.geometric(p=probex)*taurep #geometric distribution gives us number of pulses it took to get excited
            else:
                nextem = nextem + scipy.random.exponential(k_excitation)
            
            
            while nextcrash < nextem: # ignore crashes before excitation, find the next one after excitation
                if nextem - nextcrash > 5*crashtime: #if it's a really long way between, skip ahead.
                    nextcrash = nextcrash + int((nextem - nextcrash)/crashtime)*crashtime
                nextcrash = nextcrash + scipy.random.exponential(crashtime)
                while numpy.random.rand() > crashtransferprob:
                    nextcrash = nextcrash + scipy.random.exponential(crashtime)
            
            df = nextem + scipy.random.exponential(dftime)    
            #print(4)
            if numpy.random.rand() < bfrate: #this emitter had bright fission occur on this excitation
                bfem = nextem + newphoton(antibunch, k_emission) 
                    
                if numpy.random.rand() < sensitivity and bfem < min(df, nextcrash):
                        nextphoton = insert(nextphoton, bfem)
                        countbright = countbright + 1

            #regular emission
            nextem = nextem + newphoton(antibunch, k_emission)
             
            if numpy.random.rand() < sensitivity and nextem < min(df, nextcrash): # if we didn't miss it due to sensitivity
                nextphoton = insert(nextphoton, nextem)
            nextem = max(bfem, nextem)
            if nextem>df and df < nextcrash:
                nextem = nextem + max(scipy.random.exponential(tripem), scipy.random.exponential(tripem))
                # first approximation: we assume that our emitter is an dimer/aggregate with
                # multiple emitters - when it does dark fission, it becomes dark until the both of
                # these relax (i.e. cannot be excited if it is in the excited state, which we could
                # rationalize physically as the absorption changes when the molecule is in the 
                # excited state so that it cannot absorb the laser wavelength until it returns to
                # the ground state) 

        lastcrash[i] = nextcrash
        lastphoton[i] = nextem
        
        

    if numEms == 0:
        numEms = 1
        while nextOut<nextIn:
            nextOut = nextOut + d.diffuse(diffouttime, 1)
            ndiffsOut = ndiffsOut + 1
        lastphoton = [nextIn]
        lastcrash = [nextIn]
        nextIn = nextIn + d.diffuse(diffouttime, AvgEms)
        diffsIn = diffsIn + 1

    return nextphoton, lastphoton, numEms, nextIn, nextOut, diffsIn, diffsOut, ndiffsOut, countbright, lastcrash,

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
