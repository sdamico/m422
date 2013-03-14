"""
- mdct.py -- Computes reasonably fast MDCT/IMDCT using numpy FFT/IFFT
"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np

### Problem 1.a ###
def MDCTslow(data, a, b, isInverse=False):
    """
    Slow MDCT algorithm for window length a+b following pp. 130 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    (and where 2/N factor is included in forward transform instead of inverse)
    a: left half-window length
    b: right half-window length
    data is numpy array?
    """
    N=a+b
    no=float(b+1)/2
    inarr=np.array(data)
    ### YOUR CODE STARTS HERE ###
    if (isInverse):
        narr=np.array(range(N))
        karr=np.array(range(N/2))
        xn=np.array([(lambda x: np.sum(inarr*(np.cos(2*np.pi/N*(x+no)*(karr+0.5)))))(n) for n in narr])
        final=xn*float(2)
    else:
        narr=np.array(range(N))
        karr=np.array(range(N/2))
        Xk= np.array([(lambda x: np.sum(inarr*(np.cos(2*np.pi/N*(narr+no)*(x+0.5)))))(w) for w in karr])
        Xk=Xk*float(2)/N
        final=Xk
        # do a list comprehension 
    #print final
    #print final
    return final # CHANGE THIS
        ### YOUR CODE ENDS HERE ###

### Problem 1.c ###
def MDCT(data, a, b, isInverse=False):
    """
    Fast MDCT algorithm for window length a+b following pp. 141-143 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    (and where 2/N factor is included in forward transform instead of inverse)
    a: left half-window length
    b: right half-window length
    """
    
    ### YOUR CODE STARTS HERE ###

    
    inarr=np.array(data)
    N=a+b
    no=float(b+1)/2
    narr=np.array(range(N))
 
    if (isInverse):
        #need to multiply by N
        karr=np.array(range(N))
        extXk=np.concatenate((inarr, -inarr[::-1]))
        pretwid=extXk*np.exp(1j*2*np.pi/N*karr*no)
        ifftres=np.fft.ifft(pretwid)
        posttwid=np.exp(1j*2*np.pi/(2*N)*(narr+no))*ifftres
        final=np.real(posttwid)
        final=final*N #to make up for it
    else:
        karr=np.array(range(N/2))
        pretwid=np.exp(-1j*2*np.pi*narr/(2*N))
        expr1=inarr*pretwid
        fftres=np.fft.fft(expr1)
        kelem=fftres[0:N/2]
        posttwid=np.exp(-1j*2*np.pi/N*no*(karr+0.5))*kelem
        final=np.real(posttwid)
        final=final*float(2)/N   
        #return np.zeros_like( data ) # CHANGE THIS
    return final
    ### YOUR CODE ENDS HERE ###

def IMDCT(data,a,b):
    
    ### YOUR CODE STARTS HERE ###
    
    return MDCT(data,a,b,True) # CHANGE THIS
    ### YOUR CODE ENDS HERE ###

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":
    
    ### YOUR TESTING CODE STARTS HERE ###
    
    pass # THIS DOES NOTHING
    
    ### YOUR TESTING CODE ENDS HERE ###
    
