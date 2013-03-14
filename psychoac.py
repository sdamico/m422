scaleA = 100.


import numpy as np
from window import *
import pylab as py
import spectralpeaks as spk
import mdct

def SPL(intensity):
    """
    Returns the SPL corresponding to intensity (in units where 1 implies 96dB)
    """
    temp=intensity
    #intensity=np.array([(lambda x: float(x))(ii) for ii in temp])
    SPLfin=96+10*np.log(intensity)/np.log(10)
    return SPLfin
    #return np.zeros_like(intensity) # TO REPLACE WITH YOUR CODE

def Intensity(spl):
    """
    Returns the intensity (in units of the reference intensity level) for SPL spl
    """
    #temp=spl
    #spl=np.array([(lambda x: float(x))(ii) for ii in temp])
    intfin=10**((spl-96)/10)
    return intfin

def Thresh(f):
    """Returns the threshold in quiet measured in SPL at frequency f (in Hz)"""
    #temp=f
    #freq=np.array([(lambda x: float(x))(ii) for ii in temp])
    temp=f
    freq=temp/1000. #convert to kHz
    thr=3.64*(freq**(-0.8))-6.5*np.exp(-0.6*(freq-3.3)**2)+10**(-3)*(freq**4)
    return thr # TO REPLACE WITH YOUR CODE

def Bark(f):
    """Returns the bark-scale frequency for input frequency f (in Hz) """
    #temp=f
    #freq=np.array([(lambda x: float(x))(ii) for ii in temp])
    #temp=freq
    temp=f
    freq=temp/1000. #convert to kHz
    barkres=13.*np.arctan(0.76*freq)+3.5*np.arctan((freq/7.5)**2)
    return barkres # TO REPLACE WITH YOUR CODE

class Masker:
    """
    a Masker whose masking curve drops linearly in Bark beyond 0.5 Bark from the
    masker frequency
    """

    def __init__(self,f,SPL,isTonal=True):
        """
        initialized with the frequency and SPL of a masker and whether or not
        it is Tonal
        """
        self.frequency=f
        self.tonal=isTonal
        self.mspl=SPL
        #pass # TO REPLACE WITH YOUR CODE

    def IntensityAtFreq(self,freq):
        """The intensity of this masker at frequency freq"""
        z=Bark(freq)
        finalintensity=IntensityAtBark(self,z)
        return finalintensity # TO REPLACE WITH YOUR CODE

    def IntensityAtBark(self,z):
        """The intensity of this masker at Bark location z"""
        dz=z-Bark(self.frequency)
        if (self.tonal):
            delta=14.5+Bark(self.frequency)
        else:#noise masking
            delta=6
        spr=0

        if (np.abs(dz)<=0.5):
            spr=0
        elif (dz<-0.5):
            spr=-27.*(np.abs(dz)-0.5)
        else:
            spr=(-27.+0.367*np.max([self.mspl-40,0]))*(np.abs(dz)-0.5)

        finalspl=self.mspl-delta+spr
        finalintensity=Intensity(finalspl)
        return finalintensity # TO REPLACE WITH YOUR CODE

    def vIntensityAtBark(self,zVec):
        """The intensity of this masker at vector of Bark locations zVec"""
        dz=zVec-Bark(self.frequency)
        spr=np.zeros(len(zVec))
        zrvec=np.zeros(len(zVec))

        if (self.tonal):
            #delta=26
            delta=18
        else:
            delta=6

        mask1=(np.abs(dz)<=0.5)
        spr1=np.zeros(len(zVec))
        mask2= (dz<-0.5)
        spr2=-27.*(np.abs(dz)-0.5)
        mask3=(dz>0.5)
        spr3=(-27.+0.367*np.maximum(self.mspl-40,zrvec))*(np.abs(dz)-0.5)
        spr=spr1*mask1+spr2*mask2+spr3*mask3
        finalspl=(self.mspl-delta)+spr
        finalintensity=Intensity(finalspl)
        return finalintensity # TO REPLACE WITH YOUR CODE







class SpreadEnergy:
    """
    a Masker whose masking curve drops linearly in Bark beyond 0.5 Bark from the
    masker frequency
    """

    def __init__(self,f,SPL,isTonal=True):
        """
        initialized with the frequency and SPL of a masker and whether or not
        it is Tonal
        """
        self.frequency=f
        self.tonal=isTonal
        self.mspl=SPL
        #pass # TO REPLACE WITH YOUR CODE

    def IntensityAtFreq(self,freq):
        """The intensity of this masker at frequency freq"""
        z=Bark(freq)
        finalintensity=IntensityAtBark(self,z)
        return finalintensity # TO REPLACE WITH YOUR CODE

    def IntensityAtBark(self,z):
        """The intensity of this masker at Bark location z"""
        dz=z-Bark(self.frequency)
        if (self.tonal):
            delta=0
        else:#noise masking
            delta=0
        spr=0

        if (np.abs(dz)<=0.5):
            spr=0
        elif (dz<-0.5):
            spr=-27.*(np.abs(dz)-0.5)
        else:
            spr=(-27.+0.367*np.max([self.mspl-40,0]))*(np.abs(dz)-0.5)

        finalspl=self.mspl-delta+spr
        finalintensity=Intensity(finalspl)
        return finalintensity # TO REPLACE WITH YOUR CODE

    def vIntensityAtBark(self,zVec):
        """The intensity of this masker at vector of Bark locations zVec"""
        dz=zVec-Bark(self.frequency)
        spr=np.zeros(len(zVec))
        zrvec=np.zeros(len(zVec))

        if (self.tonal):
            #delta=26
            delta=0
        else:
            delta=0

        mask1=(np.abs(dz)<=0.5)
        spr1=np.zeros(len(zVec))
        mask2= (dz<-0.5)
        spr2=-27.*(np.abs(dz)-0.5)
        mask3=(dz>0.5)
        spr3=(-27.+0.367*np.maximum(self.mspl-40,zrvec))*(np.abs(dz)-0.5)
        spr=spr1*mask1+spr2*mask2+spr3*mask3
        finalspl=(self.mspl-delta)+spr
        finalintensity=Intensity(finalspl)
        return finalintensity # TO REPLACE WITH YOUR CODE








# Default data for 25 scale factor bands based on the traditional 25 critical bands
cbFreqLimits = [100,200,300,400,510,630,770,920,1080,1270,1480,1720,2000,2320,2700,3150,3700,4400,5300,6400,7700,9500,12000,15500]  # TO REPLACE WITH THE APPROPRIATE VALUES

def AssignMDCTLinesFromFreqLimits(nMDCTLines, sampleRate, flimit = cbFreqLimits):
    """
    Assigns MDCT lines to scale factor bands for given sample rate and number
    of MDCT lines using predefined frequency band cutoffs passed as an array
    in flimit (in units of Hz). If flimit isn't passed it uses the traditional
    25 Zwicker & Fastl critical bands as scale factor bands.
    """
    initfrequencyshift=0.5*float(sampleRate)/(2*nMDCTLines)
    delta=float(sampleRate)/(2*nMDCTLines)
    currentfreq=initfrequencyshift
    numMDCTlines=np.zeros(len(cbFreqLimits)+1)
    totalbandlines=0
    for ii in range(len(flimit)):
        numbandlines=0
        upperbound=flimit[ii]
        while (currentfreq<upperbound):
            numbandlines=numbandlines+1
            totalbandlines=totalbandlines+1
            currentfreq=currentfreq+delta
        numMDCTlines[ii]=numbandlines
    numMDCTlines[len(cbFreqLimits)]=nMDCTLines-totalbandlines
    return numMDCTlines
    #return np.zeros(len(flimit)) # TO REPLACE WITH YOUR CODE

class ScaleFactorBands:
    """
    A set of scale factor bands (each of which will share a scale factor and a
    mantissa bit allocation) and associated MDCT line mappings.

    Instances know the number of bands nBands; the upper and lower limits for
    each band lowerLimit[i in range(nBands)], upperLimit[i in range(nBands)];
    and the number of lines in each band nLines[i in range(nBands)]
    """

    def __init__(self,nLines):
        """
        Assigns MDCT lines to scale factor bands based on a vector of the number
        of lines in each band
        """
        self.nBands=len(nLines)
        rr=np.cumsum(nLines)
        self.lowerLine=np.concatenate((np.array([0]),np.array(rr[0:len(rr)-1])))
        self.upperLine=rr-1
        self.nLines=np.int_(np.zeros(len(nLines)))
        for iBand in range(len(nLines)):
            self.nLines[iBand]=(self.upperLine[iBand]-self.lowerLine[iBand]+1)

def CalcSMRs(data, MDCTdata, MDCTscale, sampleRate, sfBands):
    """
    Set SMR for each critical band in sfBands.

    Arguments:
                data:       is an array of N time domain samples
                MDCTdata:   is an array of N/2 MDCT frequency lines for the data
                            in data which have been scaled up by a factor
                            of 2^MDCTscale
                MDCTscale:  is an overall scale factor for the set of MDCT
                            frequency lines
                sampleRate: is the sampling rate of the time domain samples
                sfBands:    points to information about which MDCT frequency lines
                            are in which scale factor band

    Returns:
                SMR[sfBands.nBands] is the maximum signal-to-mask ratio in each
                                    scale factor band

    Logic:
                Performs an FFT of data[N] and identifies tonal and noise maskers.
                Sums their masking curves with the hearing threshold at each MDCT
                frequency location to the calculate absolute threshold at those
                points. Then determines the maximum signal-to-mask ratio within
                each critical band and returns that result in the SMR[] array.
    """
    N=len(data)
    Fs=sampleRate
    windata=HanningWindow(data) #window with a Hanning window
    xk=np.fft.fft(windata)
    xksqr=np.square(np.abs(xk))
    aintdft=float(4)/((N**2)*0.375)*xksqr
    freqint= np.arange(0,Fs,float(Fs)/(float(N)))
    frontfreq = np.float_((np.arange(len(MDCTdata))+.5)*Fs)/np.float(len(MDCTdata)*2)
    r=spk.findpeaks(aintdft)
    fronthalfdft=aintdft[0:N/2]
    freqidx=r[0]
    freqidxfloats=freqidx
    freqidx=np.array([(lambda x: int(x))(ii) for ii in freqidx]) #converts into integers
    peakintensities=r[1]
    freqidxhalf=freqidx[0:len(freqidx)/2] #peak frequency indices
    peakstrint=peakintensities[0:len(peakintensities)/2]#peak intensity values
    spldft=SPL(peakstrint)
    setzeros=np.zeros(N/2)
    freqloca=setzeros
    actualfreq=float(Fs)*(freqidxhalf)/N+ float(Fs)/N*0.5
    #create a mask for each frequency found
    maskerarray=np.array([])
    for ii in range(len(spldft)):
        m=Masker(actualfreq[ii],spldft[ii],True)
        maskerarray=np.append(maskerarray,m)
    fin=Thresh(frontfreq)
    quietthresintens=Intensity(fin)
    zvec=Bark(frontfreq)
    sumofmasks=np.array(quietthresintens) #masks are added to threshold of quiet
    for ii in range(len(maskerarray)): #addition of all the masks
        tempmask=maskerarray[ii]
        sumofmasks=sumofmasks+np.power(tempmask.vIntensityAtBark(zvec),scaleA)
    sumofmasks = np.power(sumofmasks,1./scaleA)
    maskthres=SPL(sumofmasks)

    #####do MDCT lines correspond to half frequency away?
    MDCTdatanorm=MDCTdata/float(2**MDCTscale)
    MDCTintensity= 4*np.square(MDCTdatanorm)

    SMR=SPL(MDCTintensity)-maskthres
    llim=sfBands.lowerLine
    ulim=sfBands.upperLine
    numBands=sfBands.nBands
    SMRband=np.zeros(numBands)
    for ii in range(sfBands.nBands):
        idxlower=llim[ii]
        idxhigh=ulim[ii]
        SMRband[ii]=np.max(SMR[idxlower:idxhigh+1])

    return SMRband # TO REPLACE WITH YOUR CODE

def Maskingcurve(data, freq, sampleRate):
    """
    Set SMR for each critical band in sfBands.

    Arguments:
                data:       is an array of N time domain samples
                MDCTdata:   is an array of N/2 MDCT frequency lines for the sin windowed
                            in data which have been scaled up by a factor
                            of 2^MDCTscale
                FFTdata     is an array of N
                            frequency lines
                sampleRate: is the sampling rate of the time domain samples
                sfBands:    points to information about which MDCT frequency lines
                            are in which scale factor band

    Returns:
                final= an array showing an SPL value of the masking curve for every line
    Logic:
                Performs an FFT of data[N] and identifies tonal and noise maskers.
                Sums their masking curves with the hearing threshold at each MDCT
                frequency location to the calculate absolute threshold at those
                points. Then determines the maximum signal-to-mask ratio within
                each critical band and returns that result in the SMR[] array.
    """
    N=len(data)
    Fs=sampleRate

    # compute fft
    windata=HanningWindow(data) #window with a Hanning window
    xk=np.fft.fft(windata)
    # square fft
    xksqr=np.square(np.abs(xk))
    # compute intensity
    aintdft=float(4)/((N**2)*0.375)*xksqr
    # compute fft frequency values
    freqint= np.arange(0,Fs,float(Fs)/(float(N)))
    frontfreq=freq
    # find peak intensity indices and values
    r=spk.findpeaks(aintdft)
    fronthalfdft=aintdft[0:N/2]
    # peak intensity frequency indices
    freqidx=np.int_(r[0])
    freqidxfloats=r[0]
    freqidxhalf=freqidx[0:len(freqidx)/2]
    peakintensities=r[1]
    peakstrint=peakintensities[0:len(peakintensities)/2]#peak intensity values
    spldft=SPL(peakstrint)
    setzeros=np.zeros(N/2)
    freqloca=setzeros
    actualfreq=float(Fs)*(freqidxhalf)/N+ float(Fs)/N*0.5
    #create a mask for each frequency found
    maskerarray=np.array([])
    for ii in range(len(spldft)):
        m=Masker(actualfreq[ii],spldft[ii],True)
        maskerarray=np.append(maskerarray,m)
    fin=Thresh(frontfreq)
    quietthresintens=Intensity(fin)
    zvec=Bark(freq)
    sumofmasks=np.array(quietthresintens) #masks are added to threshold of quiet
    for ii in range(len(maskerarray)): #addition of all the masks
        tempmask=maskerarray[ii]
        sumofmasks=sumofmasks+np.power(tempmask.vIntensityAtBark(zvec),scaleA)
    sumofmasks = np.power(sumofmasks,1./scaleA)
    maskthres=SPL(sumofmasks)

    return maskthres

def SpreadSignalEnergy(data, THR, freqs, sampleRate):
    """
    Set SMR for each critical band in sfBands.

    Arguments:
                data:       is an array of N time domain samples
                MDCTdata:   is an array of N/2 MDCT frequency lines for the sin windowed
                            in data which have been scaled up by a factor
                            of 2^MDCTscale
                FFTdata     is an array of N
                            frequency lines
                sampleRate: is the sampling rate of the time domain samples
                sfBands:    points to information about which MDCT frequency lines
                            are in which scale factor band

    Returns:
                final= an array showing an SPL value of the masking curve for every line
    Logic:
                Performs an FFT of data[N] and identifies tonal and noise maskers.
                Sums their masking curves with the hearing threshold at each MDCT
                frequency location to the calculate absolute threshold at those
                points. Then determines the maximum signal-to-mask ratio within
                each critical band and returns that result in the SMR[] array.
    """
    N=len(data)
    Fs=sampleRate
    windata=HanningWindow(data) #window with a Hanning window
    xk=np.fft.fft(windata)

    xksqr=np.square(np.abs(xk))
    aintdft=float(4)/((N**2)*0.375)*xksqr
    freqint= np.arange(0,Fs,float(Fs)/(float(N)))
    frontfreq=freqint[0:N/2]
    r=spk.findpeaks(aintdft)
    fronthalfdft=aintdft[0:N/2]
    freqidx=r[0]
    freqidxfloats=freqidx
    freqidx=np.array([(lambda x: int(x))(ii) for ii in freqidx]) #converts into integers
    peakintensities=r[1]
    freqidxhalf=freqidx[0:len(freqidx)/2] #peak frequency indices

    peakstrint=peakintensities[0:len(peakintensities)/2]#peak intensity values

    spldft=SPL(peakstrint)
    setzeros=np.zeros(N/2)
    freqloca=setzeros
    actualfreq=float(Fs)*(freqidxhalf)/N+ float(Fs)/N*0.5
    #create a mask for each frequency found
    maskerarray=np.array([])
    for ii in range(len(spldft)):
        m=SpreadEnergy(actualfreq[ii],spldft[ii],True)
        maskerarray=np.append(maskerarray,m)
    fin=Thresh(freqs)
    quietthresintens=Intensity(fin)
    zvec=Bark(freqs)
    sumofmasks=np.array(quietthresintens) #masks are added to threshold of quiet
    for ii in range(len(maskerarray)): #addition of all the masks
        tempmask=maskerarray[ii]
        sumofmasks=sumofmasks+np.power(tempmask.vIntensityAtBark(zvec),scaleA)
    sumofmasks = np.power(sumofmasks,1./scaleA)
    maskthres=SPL(sumofmasks)

    return maskthres # TO REPLACE WITH YOUR CODE

def spreadsignalfactor(freq):
    x=np.array([-2.5,-2.45,-2.35,-2.2,-2.1,-2,-1.8,-1.6,-1.4,-1.15,-0.8,-0.6,-0.4,-0.2,-0.1,0,0,0,0,0,0,0,0,0,0]) #values at each bark starting from 0
    x = 10*x
    z=np.arange(0,25,1)
    w=Bark(freq)
    w=np.floor(w)
    idx=w.astype(int)
    idx[idx>=len(x)]=len(x)-1
    fin=np.zeros(len(idx))
    for ii in range(len(idx)):
        fin[ii]=x[idx[ii]]
    fin=np.interp(Bark(freq),z,x)
    return fin




#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    x=np.array([0,100,200,300,400,510,630,770,920,1080,1270,1480,1720,2000,2320,2700,3150,3700,4400,5300,6400,7700,9500,12000,15500]) # TO REPLACE WITH YOUR CODE
    """
    y=Bark(x)
    print y

    y=np.arange(20,20000,1)
    fin=Thresh(y)
    r=np.log(y)/np.log(10)
    pylab.plot(r,fin)
    pylab.show()
    #print Thresh(20000)
    #print Thresh(0.02)
    """
    N=1024
    Fs=48000
    ts=np.arange(0,N)
    a0=0.45
    a1=0.15
    a2=0.15
    a3=0.1
    a4=0.06
    a5=0.05
    xn0=a0*np.cos(2*np.pi*440*ts/Fs)
    xn1=a1*np.cos(2*np.pi*550*ts/Fs)
    xn2=a2*np.cos(2*np.pi*660*ts/Fs)
    xn3=a3*np.cos(2*np.pi*880*ts/Fs)
    xn4=a4*np.cos(2*np.pi*4400*ts/Fs)
    xn5=a5*np.cos(2*np.pi*8800*ts/Fs)

    xn=xn0+xn1+xn2+xn3+xn4+xn5
    xwin=SineWindow(xn)
    xmdct=mdct.MDCT(xwin,N/2,N/2)
    linesmdct=AssignMDCTLinesFromFreqLimits(512,48000)
    sfBands=ScaleFactorBands(linesmdct)
    x=CalcSMRs(xn,xmdct,0,48000,sfBands)

    aintdft=spk.intensityDFT(xwin)
    freqint=np.arange(0,Fs,float(Fs)/(float(N)))
    frontfreq=freqint[0:N/2]
    r=findpeaks(aintdft) #find peaks
    fronthalf=aintdft[0:N/2]#only need front N/2 frequencies in DFT
    dftsig=psy.SPL(fronthalf)

    freqidx=r[0]
    freqidx=[(lambda x: int(x))(ii) for ii in freqidx]
    freqidxhalf=freqidx[0:len(freqidx)/2]
    spldft=psy.SPL(r[1])
    spldfthalf=spldft[0:len(spldft)/2]
    setzeros=np.zeros(N/2)
    freqloca=setzeros
    SPLpeaks=setzeros
    SPLpeaks[freqidxhalf]=spldfthalf
    #actualfreqpeaks=(freqlocations*float(Fs))/N
    pylab.plot(np.log(frontfreq), dftsig)
    pylab.show()

    pylab.plot(psy.Bark(frontfreq),dftsig,color='red')
    pylab.plot(psy.Bark(frontfreq),SPLpeaks,color='blue')
    zvec=psy.Bark(frontfreq)
    freq1=freqidx[0]*float(Fs)/N
    freq2=freqidx[1]*float(Fs)/N
    freq3=freqidx[2]*float(Fs)/N
    freq4=freqidx[3]*float(Fs)/N
    freq5=freqidx[4]*float(Fs)/N
    freq6=freqidx[5]*float(Fs)/N
    spl1=spldft[0]
    spl2=spldft[1]
    spl3=spldft[2]
    spl4=spldft[3]
    spl5=spldft[4]
    spl6=spldft[5]
    m1=psy.Masker(freq1,spl1,True)
    testlist=np.array([])
    testlist=np.append(testlist,m1)
    m2=psy.Masker(freq2,spl2,True)
    m3=psy.Masker(freq3,spl3,True)
    m4=psy.Masker(freq4,spl4,True)
    m5=psy.Masker(freq5,spl5,True)
    m6=psy.Masker(freq6,spl6,True)
    mask1=m1.vIntensityAtBark(zvec)
    mask2=m2.vIntensityAtBark(zvec)
    mask3=m3.vIntensityAtBark(zvec)
    mask4=m4.vIntensityAtBark(zvec)
    mask5=m5.vIntensityAtBark(zvec)
    mask6=m6.vIntensityAtBark(zvec)

    pylab.plot(zvec,psy.SPL(mask1),color='yellow')
    pylab.plot(zvec,psy.SPL(mask2),color='orange')
    pylab.plot(zvec,psy.SPL(mask3),color='green')
    pylab.plot(zvec,psy.SPL(mask4),color='pink')
    pylab.plot(zvec,psy.SPL(mask5),color='purple')
    pylab.plot(zvec,psy.SPL(mask6),color='tan')

    ######threshold of quiet
    #y=np.arange(20,20000,1)
    fin=psy.Thresh(frontfreq)
    quietthresint=psy.Intensity(fin)
    #x=psy.Bark(y)
    #r=np.log(y)/np.log(10)
    #pylab.plot(x,fin)
    #pylab.show()

    final=mask1+mask2+mask3+mask4+mask5+mask6+quietthresint
    pylab.plot(zvec,psy.SPL(final),color='black')
    #pylab.plot(zvec,fin,color='black')
    pylab.ylim([-10,100])


