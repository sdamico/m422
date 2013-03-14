"""
quantize.py -- routines to quantize and dequantize floating point values between
-1.0 and 1.0 (" signed fractions ")
"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np
import math

### Problem 1.a.i ###
def QuantizeUniform(aNum,nBits):
    """
    Uniformly quantize signed fraction aNum with nBits
    """
    #Notes:
    #The overload level of the quantizer should be 1.0

    #aQuantizedNum = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE
	### YOUR CODE STARTS HERE ###
    s=1 #0 if it is positive and 1 if it is negative
    if (aNum>=0):
        s=0
    if (abs(aNum)>=1):
        aQuantizedNum=2**(nBits-1)-1
    else:
        aQuantizedNum=math.floor((float((2**nBits-1)*abs(aNum))+1)/2)
    if (s==1):
        aQuantizedNum=int(aQuantizedNum)|(2**(nBits-1))
    else:
        aQuantizedNum=int(aQuantizedNum)
    ### YOUR CODE ENDS HERE ###

    return aQuantizedNum

### Problem 1.a.i ###
def DequantizeUniform(aQuantizedNum,nBits):
    """
    Uniformly dequantizes nBits-long number aQuantizedNum into a signed fraction
    """
    #aNum = 0.0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    s=0 #pos by defeault
    if (aQuantizedNum>=(2**(nBits-1))):
        s=1
        aQuantizedNum=aQuantizedNum-2**(nBits-1)

    aNum=(2*abs(float(aQuantizedNum)))/(2**nBits-1)
    ### YOUR CODE ENDS HERE ###
    if (s==1):
        aNum=aNum*(-1)
    return aNum

### Problem 1.a.ii ###
def vQuantizeUniform(aNumVec, nBits):
    """
    Uniformly quantize vector aNumberVec of signed fractions with nBits. Assume that the vector elements are already floats
    """
   # aQuantizedNumVec = np.zeros_like(aNumVec, dtype = int) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    #Notes:
    #Make sure to vectorize properly your function as specified in the homework instructions

    ### YOUR CODE STARTS HERE ###
    templen=np.size(aNumVec)
    afVec=np.ones(templen)
    negativeval=(aNumVec<0)
    negativebits=negativeval*(2**(nBits-1)) #will add to final array later
    temp=aNumVec
    aNumVec=np.abs(temp)#get rid of the negative signs. nNumVec now pos
    boolgreq2one=(aNumVec>=1)
    boolleq2one=(aNumVec<1)
    cutoffpart=(2**(nBits-1)-1)*boolgreq2one #will add this to final array
    noncutoff=boolleq2one*np.floor(((2**nBits-1)*np.abs(aNumVec)+1)/2)
    ### YOUR CODE ENDS HERE ###
    temp1=np.add(negativebits,cutoffpart)
    aQuantizedNumVec=np.add(temp1,noncutoff)
    tempy=aQuantizedNumVec.astype(int)
    aQuantizedNumVec=tempy #convert to integer array
    return aQuantizedNumVec


### Problem 1.a.ii ###
def vDequantizeUniform(aQuantizedNumVec, nBits):
    """
    Uniformly dequantizes vector of nBits-long numbers aQuantizedNumVec into vector of  signed fractions
    """

    #aNumVec = np.zeros_like(aQuantizedNumVec, dtype = float) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    xtemp=np.array(aQuantizedNumVec) #convert to float for division later
    ytemp=xtemp.astype(float)
    aQuantizedNumVec=ytemp
    negposarray=(aQuantizedNumVec>=(2**(nBits-1)))
    posarray=(aQuantizedNumVec<(2**(nBits-1)))
    templen=np.size(aQuantizedNumVec)
    afVec=np.ones(templen)
    afVec=afVec*aQuantizedNumVec #converts the vector into float 64 bit if it isn't yet
    afVec=afVec-(2**(nBits-1))*negposarray#converst array to positive values for now
    aNum=(2*afVec)/(2**nBits-1)
    ### YOUR CODE ENDS HERE ##
    aNumVec=aNum*posarray+(-1)*aNum*negposarray #adds the two, but uses postive and negative masks to add up correctly
    ### YOUR CODE ENDS HERE ###

    return aNumVec


###Scalefactor with the bits still. Remove later

def ScaleFactorBin(aNum, nScaleBits=3, nMantBits=5):
    """
    Return the floating-point scale factor for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """
    #Notes:
    #The scale factor should be the number of leading zeros
    #scale = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    r=2**nScaleBits-1+nMantBits-1#number of bits excluding the sign bit
    x=abs(aNum)
    binstr=bin(x)[2:].zfill(r)
    count=0#counts zeros
    for ii in binstr:
        if ii=='0':
            count=count+1
        else:
            break
    if count<(2**nScaleBits-1):
        scale=count
    else:
        scale=2**nScaleBits-1

    ### YOUR CODE ENDS HERE ###

    return scale


### Problem 1.b ###
def ScaleFactor(aNum, nScaleBits=3, nMantBits=5):
    """
    Return the floating-point scale factor for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """
    #Notes:
    #The scale factor should be the number of leading zeros
    #scale = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    r=2**nScaleBits-1+nMantBits-1#number of bits excluding the sign bit
    rtot=2**nScaleBits-1+nMantBits
    aNum=abs(aNum) #for the Scale factor, we don't care about sign so just get rid of it
    temp=QuantizeUniform(aNum,rtot)#quantize using rtot bits
    aNum=temp#aNum is now in terms of bits, where leftmost bit determines whether negative or not
    r=int(r)
    binstr=bin(aNum)[2:].zfill(r)#bitstring excluding the sign bit
    count=0#counts zeros in front
    count=len(binstr.split('1',1)[0])#counts number of zeros in front
    if count<(2**nScaleBits-1):
        scale=count
    else:
        scale=2**nScaleBits-1

    ### YOUR CODE ENDS HERE ###

    return scale

### Problem 1.b ###
def MantissaFP(aNum, scale, nScaleBits=3, nMantBits=5):
    """
    Return the floating-point mantissa for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

   # mantissa = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    if aNum==0:
        mantissa=0
    else:#not zero
        s=0 #sign
        x=abs(aNum)
        if (aNum<0):
            s=1
        r=2**nScaleBits-1+nMantBits-1
        rtot=2**nScaleBits-1+nMantBits
        quanx=QuantizeUniform(x,rtot)
        binstr=bin(quanx)[2:].zfill(r)#our binary string excluding sign bit
        if (scale==2**nScaleBits-1):
            jj=scale #set starting index to index right after scale bits
        else: #scale not equal 2**Rs-1
            jj=binstr.index('1')##index in array of first nonzero element
            jj=jj+1#actually you do not want first nonzero element
        mtemp=binstr[jj:]
        mlen=len(mtemp)
        if (mlen<nMantBits-1):
            for jj in range(nMantBits-1-mlen):
                mtemp=mtemp+'0'
        final=mtemp[0:nMantBits-1]
        mantissa=int(final,2)
        if (s==1):
            mantissa=mantissa+2**(nMantBits-1)
    ### YOUR CODE ENDS HERE ###

    return mantissa

### Problem 1.b ###
def DequantizeFP(scale, mantissa, nScaleBits=3, nMantBits=5):
    """
    Returns a  signed fraction for floating-point scale and mantissa given specified scale and mantissa bits
    """

    #aNum = 0.0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE
    s=1#when in doubt, positive
    rtot=2**nScaleBits-1+nMantBits
    ### YOUR CODE STARTS HERE ###
    if (mantissa>=2**(nMantBits-1)):#if it is negative, remove the negative bit and then store the fact that it is negative for later
        s=-1
        mantissa=mantissa-2**(nMantBits-1)
    startstr=''
    zerostr=startstr.zfill(scale)
    if (scale==(2**nScaleBits-1)):
        aNumbin=zerostr+(bin(abs(mantissa))[2:]).zfill(nMantBits-1)
    else:
        aNumbin=zerostr+'1'+(bin(abs(mantissa))[2:]).zfill(nMantBits-1)
        temp=aNumbin
        aNumbin=temp+(''.zfill(2**nScaleBits-1+nMantBits-1-len(temp)))#must append the additional zeros to the end, otherwise, you will get index out of bounds ie number 66
    ### YOUR CODE ENDS HERE ###
    if (scale<(2**nScaleBits-1-1)):#the extra minus 1 is to take into account if you have a 1 (nManbits-1) away from the right, you still don't need to add 1
        aNumbinarr=list(aNumbin)
        aNumbinarr[scale+nMantBits]='1'
        aNumbin=''.join(aNumbinarr)
    binint=int(aNumbin,2)

    aNum=DequantizeUniform(binint,rtot)
    temp=aNum*s#get the sign back in
    aNum=temp
    return aNum

### Problem 1.c.i ###
def Mantissa(aNum, scale, nScaleBits=3, nMantBits=5):
    """
    Return the block floating-point mantissa for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    #mantissa = 0 # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    if aNum==0:
        mantissa=0
    else:
        s=1
        x=abs(aNum)
        if(aNum<0):
            s=-1
        r=2**nScaleBits-1+nMantBits-1
        rtot=2**nScaleBits-1+nMantBits
        quanx=QuantizeUniform(x,rtot)
        binstr=bin(quanx)[2:].zfill(r)#our binary string, excluding sign bit
        jj=scale
        mtemp=binstr[jj:jj+nMantBits-1]
        mint=int(mtemp,2)
        if (s==-1):
            mantissa=mint| (2**(nMantBits-1)) # add neg sign bit
        else:
            mantissa=mint
    ### YOUR CODE ENDS HERE ###

    return mantissa

### Problem 1.c.i ###
def Dequantize(scale, mantissa, nScaleBits=3, nMantBits=5):
    """
    Returns a  signed fraction for block floating-point scale and mantissa given specified scale and mantissa bits
    """

    ### YOUR CODE STARTS HERE ###
    s=1 #postivie default
    if (mantissa>=2**(nMantBits-1)):
        s=-1 #negative
        mantissa=mantissa-2**(nMantBits-1)#mantissa at this point should be pos
    rtot=2**nScaleBits-1+nMantBits
    zerostr=''.zfill(scale)
    aNumbin=zerostr+(bin(mantissa)[2:]).zfill(nMantBits-1)#front values with zero filled
    temp=aNumbin
    aNumbin=temp+(''.zfill(rtot-1-len(temp)))
    if (scale+(nMantBits-1)<(rtot-1)):#made additional change rtot-1. must be <11
        aNumbinarr=list(aNumbin)
        if (not (mantissa==0)):
            aNumbinarr[scale+nMantBits-1]='1'
        aNumbin="".join(aNumbinarr)
    aNum=int(aNumbin,2)
    ### YOUR CODE ENDS HERE ###
    val=DequantizeUniform(aNum,rtot)
    temp=val*s
    aNum=temp
    return aNum

### Problem 1.c.ii ###
def vMantissa(aNumVec, scale, nScaleBits=3, nMantBits=5):
    """
    Return a vector of block floating-point mantissas for a vector of  signed fractions aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    #mantissaVec = np.zeros_like(aNumVec, dtype = int) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE

    ### YOUR CODE STARTS HERE ###
    rtot=2**nScaleBits-1+nMantBits
    binvec=vQuantizeUniform(aNumVec, rtot)
    mantvec=np.zeros(len(binvec)) #mantissa
    for ii in range(len(binvec)):
        if (binvec[ii]==0):
            mantissa=0
            mantvec[ii]=0
        else:
            s=1
            quantbits=binvec[ii]
            rtot=2**nScaleBits-1+nMantBits
            if (quantbits>=(2**(rtot-1))): #means its negative
                s=-1
                quantbits=quantbits-(2**(rtot-1))
            r=2**nScaleBits-1+nMantBits-1

            binstr=bin(int(quantbits))[2:].zfill(int(r))#our binary string, excluding sign bit
            jj=scale
            nMantBits=int(nMantBits)
            jj=int(jj)
            mtemp=binstr[jj:jj+nMantBits-1]
            if (mtemp==''):
                mtemp='0'
            mint=int(mtemp,2)
            if (s==-1):
                mantissa=mint| int(2**(int(nMantBits)-1)) # add neg sign bit
            else:
                mantissa=mint
            mantvec[ii]=mantissa
    ### YOUR CODE ENDS HERE ###

    mantissaVec=mantvec.view('int64')
    mantissaVec[:]=mantvec
    return mantissaVec

### Problem 1.c.ii ###
def vDequantize(scale, mantissaVec, nScaleBits=3, nMantBits=5):
    """
    Returns a vector of  signed fractions for block floating-point scale and vector of block floating-point mantissas given specified scale and mantissa bits
    """
    rtot=2**nScaleBits-1+nMantBits
    #aNumVec = np.zeros_like(mantissaVec, dtype = float) # REMOVE THIS LINE WHEN YOUR FUNCTION IS DONE
    valarray = np.zeros(len(mantissaVec),int) #creates integer array
    ### YOUR CODE STARTS HERE ###
    for ii in range(len(mantissaVec)):
        s=1 #postivie default
        mantissa=mantissaVec[ii]
        if (mantissa>=2**(nMantBits-1)):
            s=-1 #negative
            mantissa=mantissa-2**(nMantBits-1)#mantissa at this point should be pos
        rtot=2**nScaleBits-1+nMantBits
        scale=int(scale)
        zerostr=''.zfill(scale)
        mantissa=int(mantissa)
        aNumbin=zerostr+(bin(mantissa)[2:]).zfill(int(nMantBits)-1)#front values with zero filled
        temp=aNumbin
        aNumbin=temp+(''.zfill(int(rtot)-1-len(temp)))
        if (scale+(nMantBits-1)<(rtot-1)):#made additional change rtot-1. must be <11
            aNumbinarr=list(aNumbin)
            #if (not((mantissa==0)& scale< (2**nScaleBits-1))):
            if (not(mantissa==0)):
                #if (not(scale<(2**nScaleBits-1))):
                aNumbinarr[int(scale)+int(nMantBits)-1]='1'
            aNumbin="".join(aNumbinarr)
        aNum=int(aNumbin,2)
        if (s==-1): #negative so add in special bits. We will need for vector dequantization function later
            aNum=aNum+ (2**(rtot-1))
    ### YOUR CODE ENDS HERE ###
        #print aNum
        valarray[ii]=int(aNum)
    finalconvert=vDequantizeUniform(valarray,rtot)
    aNumVec=finalconvert
    ### YOUR CODE ENDS HERE ###

    return aNumVec

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###

    pass # THIS DOES NOTHING

    ### YOUR TESTING CODE ENDS HERE ###

