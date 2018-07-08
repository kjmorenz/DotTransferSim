import scipy
import numpy
import numpy as np

import diffuse as d

from rust_fastforward import fastforward, deltransgttem_helper

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

def alltransgttemtime(transtime, temtime):
    for i in range(len(transtime)):
        if transtime[i][0] < temtime[i][0] or transtime[i][1] < temtime[i][1]:
            return 0
    return 1

def deltransgttem(transtime, temtime):
    new_len = deltransgttem_helper(transtime, temtime)
    return transtime[:new_len], temtime[:new_len]

def allsemltallfiss(fisstime,semtime):
    for i in range(len(fisstime)):
        if fisstime[i] < semtime[i]:
            return 0
    return 1

def findmin(array):
    index_flat = numpy.argmin(array)
    index = numpy.unravel_index(index_flat, array.shape)
    return index, array[index]

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

def make_transtime_temtime(k_lexcitation, problex, nextem, taurep, nligands, k_fiss, k_sem, k_trans, k_tem, excite, probfiss):
    # Temporary variables
    nextlex = np.array([excite(k_lexcitation, problex, nextem, taurep) for _ in range(nligands)])
    nextfisstime = nextlex + scipy.random.exponential(k_fiss, nligands)
    nextsemtime = nextlex + scipy.random.exponential(k_sem, nligands)
    for i in range(nligands):
        while nextsemtime[i] < nextfisstime[i]:
            #print("3")
            nextlex[i] = excite(k_lexcitation, problex, nextsemtime[i], taurep)
            nextfisstime[i] = nextlex[i] + scipy.random.exponential(k_fiss)
            nextsemtime[i] = nextlex[i] + scipy.random.exponential(k_sem)

    fisstime = nextfisstime
    semtime = nextsemtime
    transtime = scipy.random.exponential(k_trans, (nligands, 2))
    transtime[:,0] += fisstime
    transtime[:,1] += fisstime
    temtime = scipy.random.exponential(k_tem, (nligands, 2))
    temtime[:,0] += fisstime
    temtime[:,1] += fisstime

    while True:
        #print("4")
        greater_map = transtime > temtime
        problems = np.logical_and(greater_map[:,0], greater_map[:,1])
        idx = np.where(problems)[0]

        if len(idx) == 0:
            break

        for i in idx:
            nextlex = excite(k_lexcitation, problex, np.max(temtime[i]), taurep)
            nextfisstime = nextlex + scipy.random.exponential(k_fiss)
            nextsemtime = nextlex + scipy.random.exponential(k_sem)
            while nextsemtime < nextfisstime:
                #print("5")
                nextlex = excite(k_lexcitation, problex, nextsemtime, taurep)
                nextfisstime = nextlex + scipy.random.exponential(k_fiss)
                nextsemtime = nextlex + scipy.random.exponential(k_sem)
            fisstime[i] = (nextfisstime)
            semtime[i] = (nextsemtime)

            transtime[i] = [fisstime[i] + scipy.random.exponential(k_trans),fisstime[i] + scipy.random.exponential(k_trans)]
            temtime[i] = [fisstime[i] + scipy.random.exponential(k_tem),fisstime[i] + scipy.random.exponential(k_tem)]
    return transtime, temtime


# first attempt at non-unity probability of fission
# the problem is that either all ligands fiss or don't which is a bit silly
def make_transtime_temtime_f(k_lexcitation, problex, nextem, taurep, nligands, k_fiss, k_sem, k_trans, k_tem, excite, probfiss):
    # Temporary variables
    nextlex = np.array([excite(k_lexcitation, problex, nextem, taurep) for _ in range(nligands)])
    
    if numpy.random.rand() < probfiss:#either all fiss or all don't
        nextfisstime = nextlex + scipy.random.exponential(k_fiss, nligands)
        nextsemtime = nextlex + scipy.random.exponential(k_sem, nligands)
        for i in range(nligands):
            while nextsemtime[i] < nextfisstime[i]:
                nextlex[i] = excite(k_lexcitation, problex, nextsemtime[i], taurep)
                nextfisstime[i] = nextlex[i] + scipy.random.exponential(k_fiss)
                nextsemtime[i] = nextlex[i] + scipy.random.exponential(k_sem)

        fisstime = nextfisstime
        semtime = nextsemtime
    
        transtime = scipy.random.exponential(k_trans, (nligands, 2))
        transtime[:,0] += fisstime
        transtime[:,1] += fisstime
        temtime = scipy.random.exponential(k_tem, (nligands, 2))
        temtime[:,0] += fisstime
        temtime[:,1] += fisstime
        while True:
            greater_map = transtime > temtime
            problems = np.logical_and(greater_map[:,0], greater_map[:,1])
            idx = np.where(problems)[0]

            if len(idx) == 0:
                break

            for i in idx:
                nextlex = excite(k_lexcitation, problex, np.max(temtime[i]), taurep)
                nextfisstime = nextlex + scipy.random.exponential(k_fiss)
                nextsemtime = nextlex + scipy.random.exponential(k_sem)
                while nextsemtime < nextfisstime:
                    #print("5")
                    nextlex = excite(k_lexcitation, problex, nextsemtime, taurep)
                    nextfisstime = nextlex + scipy.random.exponential(k_fiss)
                    nextsemtime = nextlex + scipy.random.exponential(k_sem)
                fisstime[i] = (nextfisstime)
                semtime[i] = (nextsemtime)

                transtime[i] = [fisstime[i] + scipy.random.exponential(k_trans),fisstime[i] + scipy.random.exponential(k_trans)]
                temtime[i] = [fisstime[i] + scipy.random.exponential(k_tem),fisstime[i] + scipy.random.exponential(k_tem)]

    else:
        transtime = scipy.random.exponential(k_trans, (nligands, 1))
        temtime = scipy.random.exponential(k_sem, (nligands, 1)) #temtime actually holds the singlet emission time in the case of no fission
        while True:
            greater_map = transtime > temtime
            problems = np.logical_and(greater_map[:,0], greater_map[:,1])
            idx = np.where(problems)[0]

            if len(idx) == 0:
                break

            for i in idx:
                nextlex = excite(k_lexcitation, problex, np.max(temtime[i]), taurep)
                nextfisstime = nextlex + scipy.random.exponential(k_fiss)
                nextsemtime = nextlex + scipy.random.exponential(k_sem)
                while nextsemtime < nextfisstime:
                    #print("5")
                    nextlex = excite(k_lexcitation, problex, nextsemtime, taurep)
                    nextfisstime = nextlex + scipy.random.exponential(k_fiss)
                    nextsemtime = nextlex + scipy.random.exponential(k_sem)
                fisstime[i] = (nextfisstime)
                semtime[i] = (nextsemtime)

                transtime[i] = [fisstime[i] + scipy.random.exponential(k_trans),fisstime[i] + scipy.random.exponential(k_trans)]
                temtime[i] = [fisstime[i] + scipy.random.exponential(k_tem),fisstime[i] + scipy.random.exponential(k_tem)]
    return transtime, temtime

def nextphotonss(lastphoton, sensitivity, nligands,
               k_demission, k_fiss, k_trans, k_sem, k_tem, k_dexcitation, k_lexcitation,
               diffouttime, numEms, nextOut, nextIn, endround,
               antibunch, pulsed, taurep,
               probdex,problex, AvgEms, diffsIn,
               diffsOut, ndiffsOut, testdummy, nextdex, seq, probfiss):

    #Note that as it stands, this only works for:
    # 1) PULSED - could be fixed by allowing the dot to be re-excited after emitting
    # before another triplet transfer. We would also have to consider new ligand excitations.

    ##print("c")
    make = make_transtime_temtime
    if probfiss < 1:
        make = make_transtime_temtime_f
    nextphoton = []

    newphoton = newcphoton
    if pulsed == 1:
        excite = excitepulsed
    else:
        excite = excitects

    lastnextem = 0
    lastnextdex = 0
    
    while nextOut < endround: #calculate emissions before it diffuses out
        #print("1no")
        diffOut = numpy.random.randint(numEms)#index identity of emitter which diffused out
        nextem = lastphoton[diffOut] #lastphoton is an array of length numEms holding the time of last photon emitted for each emitter

        while nextem < nextOut:
            #print("2no")
            transtime, temtime = make(k_lexcitation, problex, nextem, taurep, nligands, k_fiss, k_sem, k_trans, k_tem, excite, probfiss)
            count = 0
            while transtime.shape[0] != 0:
                count = count + 1
                if count > 110:
                    print(transtime.shape[0])
                if count > 150:
                    print(transtime)
                    print(nextdex)
                    print(nextem)
                    print(minval)
                    print(nextOut)
                    print(nextphoton)
                    crash = please
                #print("6")
                '''if lastnextem == nextem and lastnextdex == nextdex:
                    #print("minval = " + str(minval))
                    #print("nextem = " + str(nextem))
                    #print("nextOut = " + str(nextOut))
                    #print("nextdex = " + str(nextdex))
                    for i in range(len(transtime)):
                        #print("transtime i = ")
                        #print(transtime[i])
                        #print("temtime i = ")
                        #print(temtime[i])'''

                index, minval = findmin(transtime)
                if nextdex < nextem:
                    #print("7")
                    #re-excite the dot and find it's emission time
                    nextdex = excite(k_dexcitation, probdex, nextem, taurep)

                while nextdex < minval: #usually nextdex is waaayyy later
                    #print("8")
                    nextem = nextdex + scipy.random.exponential(k_demission)
                    if nextdex > nextOut or nextem > nextOut:
                        transtime = None
                        break
                    if nextem > minval:
                        nextem = minval
                        transtime[index] = float("inf")
                    elif numpy.random.rand() < sensitivity and nextem < nextOut: # if we didn't miss it due to sensitivity
                        nextphoton = insert(nextphoton, nextem)
                    
                    nextdex = excite(k_dexcitation, probdex, nextem, taurep)
                    transtime, temtime = deltransgttem(transtime,temtime)
                    if transtime.shape[0] == 0:
                        transtime = None
                        break
                    index, minval = findmin(transtime)
                if minval > nextOut or transtime is None:
                    transtime = None
                    break

                if minval >= nextem: #what if minval, or all vals, < nextem?
                    ##print(8.5)
                    nextem = minval + scipy.random.exponential(k_demission)#excitation transfers to dot and emits
                    transtime[index] = float("inf")
                    transtime, temtime = deltransgttem(transtime,temtime)
                    if transtime.shape[0] == 0:
                        ##print("hit the end")
                        transtime = None
                        if numpy.random.rand() < sensitivity and nextem < nextOut: # always add last photon in list
                            nextphoton = insert(nextphoton, nextem)
                        break
                    index, minval = findmin(transtime)
                    if nextem > minval:
                        nextem = minval
                        transtime[index] = float("inf")
                        transtime, temtime = deltransgttem(transtime,temtime)
                    elif numpy.random.rand() < sensitivity and nextem < nextOut: # if we didn't miss it due to sensitivity
                        nextphoton = insert(nextphoton, nextem)
                elif minval <= nextem:#pretty sure this is impossible
                    print("minval <= nextem")
                    print(minval)
                    print(nextem)
                    print(transtime)
                    print(nextdex)
                    crash = please
            lastnextem = nextem
            lastnextdex = nextdex
            ##print(nextem)
            ##print(nextOut)
            ##print(nextdex)


        numEms = numEms - 1
        del lastphoton[diffOut]
        diffsOut = diffsOut + 1
        nextOut = nextOut + d.diffuse(diffouttime, numEms)
        if numEms == 0:
            break #we'll want to start a new round

    while nextIn < endround:
        #print("9")
        lastphoton.append(nextIn)
        diffsIn = diffsIn + 1
        numEms = numEms + 1
        nextIn = nextIn + d.diffuse(diffouttime, AvgEms)
        if numEms == 1:
            endround = nextIn
            break #skip to next round if we just diffused out and back in'''




    for j in range(numEms):
        nextem = lastphoton[j]
        ##print("10")
        
        while nextem < endround: #continue adding photons until end of round, on average 10 emissions
            #print("11")
            transtime, temtime = make(k_lexcitation, problex, nextem, taurep, nligands, k_fiss, k_sem, k_trans, k_tem, excite, probfiss)
            count = 0
            while transtime.shape[0] != 0:
                count = count + 1
                if count > 110:
                    print(transtime.shape[0])
                if count > 150:
                    print(transtime)
                if count > 200:
                    print(nextdex)
                    print(nextem)
                    print(minval)
                    print(nextOut)
                    print(nextphoton)
                    crash = please
                #print("12")
                ##print("transtime isn't empty")
                ##print(transtime)
                index, minval = findmin(transtime)
                if nextdex < nextem:
                    #print("13")
                    #re-excite the dot and find it's emission time
                    nextdex = excite(k_dexcitation, probdex, nextem, taurep)

                while nextdex < minval: #usually nextdex is waaayyy later
                    #print("14")
                    nextem = nextdex + scipy.random.exponential(k_demission)
                    if nextem > minval:
                        nextem = minval
                        transtime[index] = float("inf")
                    elif numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                        nextphoton = insert(nextphoton, nextem)
                        ##print("3")
                    nextdex = excite(k_dexcitation, probdex, nextem, taurep)
                    transtime, temtime = deltransgttem(transtime,temtime)
                    if transtime.shape[0] == 0:
                        transtime = None
                        break
                    index, minval = findmin(transtime)
                
                if transtime is None or transtime.shape[0] == 0:
                    #transtime = None
                    break

                if minval >= nextem:#what if minval < nextem
                    ##print(14.5)
                    nextem = minval + scipy.random.exponential(k_demission)#excitation transfers to dot and emits
                    transtime[index] = float("inf")
                    transtime, temtime = deltransgttem(transtime,temtime)
                    if transtime.shape[0] == 0:
                        ##print("otherhittheend")
                        transtime = None
                        if numpy.random.rand() < sensitivity and nextem < nextOut: # always add last photon in list
                            nextphoton = insert(nextphoton, nextem)
                        break
                    index, minval = findmin(transtime)
                    if nextem > minval:
                        nextem = minval
                        transtime[index] = float("inf")
                        transtime, temtime = deltransgttem(transtime,temtime)
                    elif numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                        nextphoton = insert(nextphoton, nextem)
                        ##print("4")
                if minval < nextem:
                    print("err")
        lastphoton[j] = nextem


    if numEms == 0:
        numEms = 1
        while nextOut<nextIn:
            ##print("15")
            nextOut = nextOut + d.diffuse(diffouttime, 1)
            ndiffsOut = ndiffsOut + 1
        lastphoton = [nextIn]
        nextIn = nextIn + d.diffuse(diffouttime, AvgEms)
        diffsIn = diffsIn + 1

    return nextphoton, lastphoton, numEms, nextIn, nextOut, diffsIn, diffsOut, ndiffsOut, nextdex

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
