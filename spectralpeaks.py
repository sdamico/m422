import numpy as np
import pylab
import window as wn
import mdct
import psychoac as psy
import matplotlib.pyplot as plt
"""
#input numpy array, assume its floats
def peaksum(xn):
    y=xn*np.conj(xn)
    print y
    temp=np.real(y)
    y=temp
    print y
    fin=np.zeros(len(y))
    for ii in range(len(y)):
        if (ii==0):
            fin[ii]=np.sum(y[ii:ii+2])
        elif (ii==(len(fin)-1)):
            fin[ii]=np.sum(y[ii-1:ii+1])
        else:
            fin[ii]=np.sum(y[ii-1:ii+2])
    return fin
"""
#input Xk squared and also index array idx acquired from mdctfindpeaksidx
def mdctfindpeaksval(xin,idx):
    val=np.array([])
    for ii in idx:
        if ((ii-1>=0) & (ii+1<len(xin))):
            temp=xin[ii]+xin[ii-1]+xin[ii+1]
            val=np.append(val,temp)
        elif (ii-1<0):
            temp=xin[ii]+xin[ii+1]
            val=np.append(val,temp)
        else:
            temp=xin[ii]+xin[ii-1]
            val=np.append(val,temp)
    return val



def mdctfindpeaksidx(xin):
    fin=np.array([])
    idx=np.array([])
    val=np.array([])
    for ii in (range(len(xin))):
        if ((ii-3>=0) & (ii+4<=(len(xin)-1))):
            if (xin[ii]==np.max(xin[ii-3:ii+4])):
                temp=xin[ii]
                #val=np.append(val,temp)
                idx=np.append(idx,ii)
        elif (ii-3<0):
            if (xin[ii]==np.max(xin[0:ii+4])):
                temp=xin[ii]
                #val=np.append(val,temp)
                idx=np.append(idx,ii)
        else:
            if (xin[ii]==np.max(xin[(ii-3):len(xin)])):
                idx=np.append(idx,ii)
    return idx

def findpeaks(xin):
    fin=np.array([])
    idx=np.array([])
    val=np.array([])
    for ii in (range(len(xin))):
        if (ii==0):
            if (xin[ii+1]-xin[ii]<0):
                temp=xin[ii]+xin[ii+1]
                val=np.append(val,temp)
                idx=np.append(idx,ii)
        elif (ii==(len(xin)-1)):
            if ((xin[ii]-xin[ii-1])>0):
                temp=xin[ii-1]+xin[ii]
                val=np.append(val,temp)
                idx=np.append(idx,ii)
        else:
            if(((xin[ii]-xin[ii-1])>0) & ((xin[ii+1]-xin[ii])<0)):
                temp=xin[ii-1]+xin[ii]+xin[ii+1]
                val=np.append(val,temp)
                idx=np.append(idx,ii)
    fin=[idx,val]
    return fin

#take in a signal, float numpy vector calculate DFT and find intensity based on that
#don't forget your signal must be windowed ahead of time

def intensityDFT(xin):
    xk=np.fft.fft(xin)
    N=len(xin)
    xksqr=np.square(np.abs(xk))
    intensefin=float(4)/((N**2)*0.375)*xksqr
    return intensefin

#take in signal, float numpy vector
#N should be even
def intensityMDCT(xin):
    N=len(xin)
    xk=mdct.MDCT(xin,N/2,N/2)
    xksqr=xk**2
    intensefin=float(8)/((N**2)*0.375)*xksqr
    return intensefin

#Testing code
if __name__ == "__main__":

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
    ####################################DFT of signal
    xn=xn0+xn1+xn2+xn3+xn4+xn5
    xnwin=wn.HanningWindow(xn)
    aintdft=intensityDFT(xnwin)
    freqint=np.arange(0,Fs,float(Fs)/(float(N)))
    frontfreq=np.array(freqint[0:N/2])
    r=findpeaks(aintdft) #find peaks
    fronthalf=aintdft[0:N/2]#only need front N/2 frequencies in DFT
    dftsig=psy.SPL(fronthalf)
    testfreq=float(Fs)*frontfreq/N

    freqidx=r[0]
    freqidx=[(lambda x: int(x))(ii) for ii in freqidx]
    freqidxhalf=np.array(freqidx[0:len(freqidx)/2])
    spldft=psy.SPL(r[1])
    spldfthalf=spldft[0:len(spldft)/2]
    setzeros=np.zeros(N/2)
    freqloca=setzeros
    SPLpeaks=setzeros
    SPLpeaks[freqidxhalf]=spldfthalf
    #actualfreqpeaks=(freqlocations*float(Fs))/N
    tempvar=frontfreq*float(Fs)/N
    pylab.plot(np.log(frontfreq)/np.log(10), dftsig)
    pylab.xlabel('log of frequency')
    pylab.ylabel('SPL')
    pylab.show()

    pylab.plot(np.log(frontfreq)/np.log(10), dftsig)
    pylab.plot(np.log(frontfreq)/np.log(10),psy.Thresh(frontfreq))
    pylab.xlabel('log of frequency')
    pylab.ylabel('SPL')

    pylab.show()
    pylab.plot(psy.Bark(frontfreq),dftsig,color='red')
    #pylab.plot(psy.Bark(frontfreq),SPLpeaks,color='blue')
    #pylab.plot(psy.Bark(frontfreq),psy.Thresh(frontfreq),color='gray')
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

    cbFreqLimits = np.array([100.,200,300,400,510,630,770,920,1080,1270,1480,1720,2000,2320,2700,3150,3700,4400,5300,6400,7700,9500,12000,15500])  # TO REPLACE WITH THE APPROPRIATE VALUES
    #for ii in cbFreqLimits:
    #    plt.axvline(x=psy.Bark(ii),color='lime')
    pylab.xlabel('bark')
    pylab.ylabel('SPL')


    pylab.ylim([-10,100])
    w=psy.Maskingcurve(xn, Fs, frontfreq)
    pylab.plot(zvec,w,color='blue')
    pylab.show()


    """ test code for peaksum
    x=np.array([1,2,3,4,5])
    y=peaksum(x)
    print x
    print y
    """
    #xk=np.fft.fft(np.abs(xk))
    #pylab.plot(xk)
    #pylab.show()

