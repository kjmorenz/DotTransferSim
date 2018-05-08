import numpy
def diffuse(diffouttime, Ems):
    nextev = float('inf')
    t = nextev
    num = round(Ems)
    if num < 1:
        num = 1 #need to fix this for very low average emitters, but generally don't have <0.5 average emitters
        
    for i in range(num):
        while t <=0 or t == float('inf'):
            t = -diffouttime*numpy.log(1-numpy.random.rand()) #poissonian with average time diffouttime
            
        if nextev > t:
            nextev = t
        t = float('inf')
        
    return nextev

def newNumEms(numEms, AvgEms, nextIn, nextOut, diffouttime, sigma): #Not used in sim but here for testing purposes
    nextdiff = min(nextIn, nextOut)
    if nextdiff == nextIn:
        if nextdiff == nextOut:
            nextIn = nextIn + d.diffuse(diffouttime, AvgEms)
            nextOut = nextOut + d.diffuse(diffouttime, numEms)

        else:
            numEms = numEms + 1
            nextIn = nextIn + diffuse(diffouttime,AvgEms)
    else:
        if numEms > 1:
            numEms = numEms - 1
            nextOut = nextOut + diffuse(diffouttime,2**numEms)
        else:
            while nextOut < nextIn:
                nextOut = nextOut + diffuse(diffouttime, 1)
            nextIn = nextIn + diffuse(diffouttime,AvgEms)


    return nextIn, nextOut, numEms
