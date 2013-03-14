from pcmfile import *
import quantize as qz
import window as wn
import mdct
import psychoac as psy
import bitalloc as ball
# Testing the PCMFile class
if __name__=="__main__":

    # create the audio file objects of the appropriate audioFile type
    #inFile= PCMFile("HW5_Test_File.wav")
    inFile= PCMFile("HW5_testFile.wav")
    outFile = PCMFile("hw5output.wav")
    Fs=48000
    # open input file and get its coding parameters
    codingParams= inFile.OpenForReading()
    cbFreqLimits = [100,200,300,400,510,630,770,920,1080,1270,1480,1720,2000,2320,2700,3150,3700,4400,5300,6400,7700,9500,12000,15500]  # TO REPLACE WITH THE APPROPRIATE VALUES

    # set additional coding parameters that are needed for encoding/decoding
    codingParams.nSamplesPerBlock = 1024

    # open the output file for writing, passing needed format/data parameters
    outFile.OpenForWriting(codingParams)

    nSampB=codingParams.nSamplesPerBlock
    nChannels=codingParams.nChannels
    # Read the input file and pass its data to the output file to be written
    ii=0
    # 1161 and 2526
    bitsleft=1161
    priorBlock=[[0 for ii in range(nSampB)] for iCh in range(nChannels)]
    overlapAndAdd=[[0 for ii in range(nSampB)] for iCh in range(nChannels)]
    numlarger=0 #number of block smaller
    firstpass=0
    starting= codingParams.numSamples
    totalnum=0
    while True:
        data=inFile.ReadDataBlock(codingParams)
        if not data: break  # we hit the end of the input file
        nBits = 8
        subsize=16
        count=0
        flag=0
        flag2=0
        for iCh in range(codingParams.nChannels):
            current=data[iCh]
            if (len(current)>1024):
                numlarger=numlarger+1
            combine=np.concatenate((priorBlock[iCh],current)) #should be concat option

            priorBlock[iCh]=current
            res=wn.SineWindow(combine)
            #print len(res)
            ######################################################## MDCT
            mdcttemp=mdct.MDCT(res,nSampB,nSampB)
            #if (flag==0):#this if statement is for debugging
            #    #print mdcttemp
            #    t1=mdcttemp

            mylines=psy.AssignMDCTLinesFromFreqLimits(1024,Fs,cbFreqLimits)
            mybands=psy.ScaleFactorBands(mylines)

            #print len(combine)
            #print len(res)
            smrmat=psy.CalcSMRs(combine,mdcttemp,0,Fs,mybands)

            #1161 or 2526
            #nmantbits=ball.BitAllocUniform(bitsleft, 10,25,mylines,smrmat) #these are the mantissa bits array
            #nmantbits=ball.BitAllocUniform(1600,10,25,mylines,smrmat)
            #nmantbits=ball.BitAllocConstSNR(bitsleft,10,25,mylines,smrmat)
            #nmantbits=ball.BitAllocConstMNR(bitsleft,16,25,mylines,smrmat)
            nmantbits=ball.BitAlloc(bitsleft, 16, 25,mylines,smrmat)
            vMant=np.zeros(len(mdcttemp)) #one perline
            vscale=np.zeros(mybands.nBands)   #one per band
    #####quantization process
            for ii in range(mybands.nBands): #loop through the bands
                rr=mdcttemp[mybands.lowerLine[ii]:mybands.upperLine[ii]+1] #get array of elements within that band

                scaleforblock=np.min(np.array([(lambda x: qz.ScaleFactor(x,4,nmantbits[ii]))(w) for w in rr])) #get the max scale within the subblockk
                #if (flag2==0): #debuggs whether the max of scalefactor is correct
                 #    t6=np.array([(lambda x: qz.ScaleFactor(x,4,nmantbits[ii]))(w) for w in rr])
                #    flag2=1

                #print scaleforblock
                vscale[ii]=scaleforblock
                vMant[mybands.lowerLine[ii]:mybands.upperLine[ii]+1]=qz.vMantissa(mdcttemp[mybands.lowerLine[ii]:mybands.upperLine[ii]+1],scaleforblock,4,nmantbits[ii])
            ######dequantization process
            #if (flag==0): #debugging whether the output scale, mantissa are good
            #    t3=vscale
            #    t4=vMant
            xfin=np.zeros(len(mdcttemp)) #xfin is our guy after dequantization
            #print xfin.dtype
            #print vMant.dtype
            for ii in range(mybands.nBands):
                xfin[mybands.lowerLine[ii]:mybands.upperLine[ii]+1]=qz.vDequantize(vscale[ii],vMant[mybands.lowerLine[ii]:mybands.upperLine[ii]+1],4,nmantbits[ii])

            #print vMant
            #if (flag==0): #debugging, displays the quantized version
            #    t2=xfin
            #    flag=1

            mdctq = xfin
            ################################################ inverse MDCT
            mdcttemp2= mdct.IMDCT(mdctq,nSampB,nSampB)

            res2=wn.SineWindow(mdcttemp2)
            leftside=res2[0:nSampB]
            rightsideov=overlapAndAdd[iCh]
            finalsum=[(leftside[ii]+rightsideov[ii]) for ii in range(len(rightsideov))]
            if (finalsum[5]>0):
                count=count+1
            overlapAndAdd[iCh]=res2[nSampB:2*nSampB]
            data[iCh] = np.array(finalsum)
            if (iCh==0):
                totalnum=totalnum+len(finalsum)
        outFile.WriteDataBlock(data,codingParams)
        #print ii
        #ii=ii+1
    for iCh in range(codingParams.nChannels): #takes care of the final block
        data=[ np.zeros(nSampB) for ii in range(codingParams.nChannels)]
        lastbuf=np.zeros(nSampB)
        combine=np.concatenate((priorBlock[iCh],lastbuf))
        res=wn.SineWindow(combine)

        mdcttemp=mdct.MDCT(res,nSampB,nSampB)

        mylines=psy.AssignMDCTLinesFromFreqLimits(1024,Fs,cbFreqLimits)
        mybands=psy.ScaleFactorBands(mylines)

        #print len(combine)
        #print len(res)
        smrmat=psy.CalcSMRs(combine,mdcttemp,0,Fs,mybands)
        #1161 or 2526
        #nmantbits=ball.BitAllocUniform(bitsleft, 10,25,mylines,smrmat) #these are the mantissa bits array
        #nmantbits=ball.BitAllocUniform(1600,10,25,mylines,smrmat)
        #nmantbits=ball.BitAllocConstSNR(bitsleft,10,25,mylines,smrmat)
        #nmantbits=ball.BitAllocConstMNR(bitsleft,16,25,mylines,smrmat)
        nmantbits=ball.BitAlloc(bitsleft, 16, 25,mylines,smrmat)
        vMant=np.zeros(len(mdcttemp)) #one perline
        vscale=np.zeros(mybands.nBands)   #one per band
        #####quantization process
        for ii in range(mybands.nBands): #loop through the bands
            rr=mdcttemp[mybands.lowerLine[ii]:mybands.upperLine[ii]+1] #get array of elements within that band

            scaleforblock=np.min(np.array([(lambda x: qz.ScaleFactor(x,4,nmantbits[ii]))(w) for w in rr])) #get the max scale within the subblockk
            #if (flag2==0): #debuggs whether the max of scalefactor is correct
            #    t6=np.array([(lambda x: qz.ScaleFactor(x,4,nmantbits[ii]))(w) for w in rr])
            #    flag2=1

            #print scaleforblock
            vscale[ii]=scaleforblock
            vMant[mybands.lowerLine[ii]:mybands.upperLine[ii]+1]=qz.vMantissa(mdcttemp[mybands.lowerLine[ii]:mybands.upperLine[ii]+1],scaleforblock,4,nmantbits[ii])
            ######dequantization process
            #if (flag==0): #debugging whether the output scale, mantissa are good
            #    t3=vscale
            #    t4=vMant
        xfin=np.zeros(len(mdcttemp)) #xfin is our guy after dequantization
            #print xfin.dtype
            #print vMant.dtype
        for ii in range(mybands.nBands):
            xfin[mybands.lowerLine[ii]:mybands.upperLine[ii]+1]=qz.vDequantize(vscale[ii],vMant[mybands.lowerLine[ii]:mybands.upperLine[ii]+1],4,nmantbits[ii])

            #print vMant
            #if (flag==0): #debugging, displays the quantized version
            #    t2=xfin
            #    flag=1

        mdctq = xfin

        mdcttemp2= mdct.IMDCT(mdctq,nSampB,nSampB)

        res2=wn.SineWindow(mdcttemp2)
        leftside=res2[0:nSampB]
        rightsideov=overlapAndAdd[iCh]
        finalsum=[(leftside[ii]+rightsideov[ii]) for ii in range(len(rightsideov))]
        if (finalsum[5]>0):
            count=count+1
        if (iCh==0):
            totalnum=totalnum+len(finalsum)
        data[iCh]=np.array(finalsum)
    #print codingParams.numSamples
    codingParams.numSamples=totalnum #final number of samples
    outFile.WriteDataBlock(data,codingParams)

    # end loop over reading/writing the blocks

    # close the files
    ending=codingParams.numSamples

    inFile.Close(codingParams)
    outFile.Close(codingParams)
    #print t1[0:20] #mdct before quantization
    #print t2[0:20] #mdct after quantization
    #print t3[0:5] #scale factors for each band
    #print t4[0:20] #mantissa for each
    #print t6 #scalefactors for that first band
    #print mylines
    #print mybands.lowerLine
    #print mybands.upperLine
