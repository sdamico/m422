"""
window.py -- Defines functions to window an array of data samples
"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np
import pylab
### Problem 1.d ###
def SineWindow(dataSampleArray):
    """
    Returns a copy of the dataSampleArray sine-windowed
    Sine window is defined following pp. 106-107 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    """
    #assumes input array are floats
    ### YOUR CODE STARTS HERE ###
    idata=np.array(dataSampleArray)
    N=len(dataSampleArray)
    x=np.array(range(N))
    window=np.sin(np.pi*(x+0.5)/N)
    final=idata*window
    return final # CHANGE THIS
    ### YOUR CODE ENDS HERE ###


def Sinedeform(dataSampleArray,isReverse=False):
    """
    Returns a copy of the dataSampleArray sine-windowed
    Sine window is defined following pp. 106-107 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    """
    Nlong=2048
    Nshort=256
    temp=np.ones(len(dataSampleArray))
    temp=temp.astype(float)
    starting=np.zeros(len(dataSampleArray))
    mask=starting
    mask[0:Nlong/2]=1
    print "this is mask"
    print "this is starting"
    firsthalf=SineWindow(temp)*mask
    #assumes input array are floats
    ### YOUR CODE STARTS HERE ###
    secondpart=np.zeros_like(starting)
    secondpart[Nlong/2:(Nlong/2+Nlong/4-Nshort/4)]=1
    thirdpart=np.zeros_like(starting)

    thirdpart[(3./4*Nlong-Nshort/4):(3./4*Nlong-Nshort/4+Nshort/2)]=1
    x=np.array(range(Nlong))
    thrsin=np.sin(np.pi*(x-(3.*Nlong/4-3*Nshort/4))/(Nshort))

    thirdpart=thirdpart*thrsin
    idata=np.array(dataSampleArray)
    N=len(dataSampleArray)
    x=np.array(range(N))
    window=np.sin(np.pi*(x+0.5)/N)
    #window=
    final=firsthalf+secondpart+thirdpart

    if (isReverse):
        return final[::-1]
    else:
        return final
    ### YOUR CODE ENDS HERE ###


def shortwind(dataSampleArray,ii=0):
    idata=np.array(dataSampleArray)
    Nshort=256
    Nlong=len(dataSampleArray)
    maskingwindow=np.zeros(Nlong)
    numelem=Nlong/Nshort
    window=np.zeros(len(dataSampleArray))
    x=np.array(range(Nlong))
    #for ii in range((len(dataSampleArray)/Nshort)):
    maskingwindow=np.zeros(Nlong) #reset dude to zero
    if (ii==0):
        maskingwindow[0:(-Nshort/4+2*Nshort/2-Nshort/4)]=1
    else:
        maskingwindow[ii*Nshort/2-Nshort/4:(ii+2)*Nshort/2-Nshort/4]=1
    window=window+np.sin(np.pi*(x-Nshort/2*ii+0.5+Nshort/4)/(Nshort))*maskingwindow
    window=window*dataSampleArray
    return window

def HanningWindow(dataSampleArray):
    """
    Returns a copy of the dataSampleArray Hanning-windowed
    Hann window is defined following pp. 106-107 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    """

    ### YOUR CODE STARTS HERE ###
    idata=np.array(dataSampleArray)
    N=len(dataSampleArray)
    x=np.array(range(N))
    window=0.5*(1-np.cos(2*np.pi*(x+0.5)/N))
    final=idata*window
    return final # CHANGE THIS
    ### YOUR CODE ENDS HERE ###



### Problem 1.d - OPTIONAL ###
def KBDWindow(dataSampleArray,alpha=4.):
    """
    Returns a copy of the dataSampleArray KBD-windowed
    KBD window is defined following pp. 108-109 and pp. 117-118 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    """

    ### YOUR CODE STARTS HERE ###

    return np.zeros_like(dataSampleArray) # CHANGE THIS
    ### YOUR CODE ENDS HERE ###

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###

    pass # THIS DOES NOTHING

    ### YOUR TESTING CODE ENDS HERE ###

