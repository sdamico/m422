"""
codec.py -- The actual encode/decode functions for the perceptual audio codec

-----------------------------------------------------------------------
© 2009 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------
2013: Edited by Sam D'Amico, Wendy Nie, Jonathan Lu
"""

import pylab as py

import numpy as np  # used for arrays

# used by Encode and Decode
import window as wn
from mdct import MDCT,IMDCT  # fast MDCT implementation (uses numpy FFT)
from quantize import *  # using vectorized versions (to use normal versions, uncomment lines 18,67 below defining vMantissa and vDequantize)

# used only by Encode
import psychoac as psy # calculates SMRs for each scale factor band
from bitalloc import *  #allocates bits to scale factor bands given SMRs
import matplotlib.pyplot as plt

def DecodeStereo(scaleFactor,bitAlloc,mantissa,overallScaleFactor,codingParams):
    """Reconstitutes a single-channel block of encoded data into a block of
    signed-fraction data based on the parameters in a PACFile object"""

    rescaleLevel = 1.*float(1<<overallScaleFactor)
    halfN = codingParams.nMDCTLines
    N = 2*halfN

    mylines=psy.AssignMDCTLinesFromFreqLimits(codingParams.nMDCTLines, codingParams.sampleRate)

    # reconstitute the first halfN MDCT lines of this channel from the stored data
    PreMDCTLine = np.zeros((2,halfN),dtype=np.float64)
    iMant=0
    # extract coded values
    for ii in range(2):
        iMant=0
        for iBand in range(codingParams.sfBands.nBands):
            nLines=mylines[iBand]
            if bitAlloc[ii][iBand]:
                PreMDCTLine[ii][iMant:(iMant+nLines)]=vDequantize(scaleFactor[ii][iBand], mantissa[ii][iMant:(iMant+nLines)],codingParams.nScaleBits, bitAlloc[ii][iBand])
            iMant += nLines
    lr=codingParams.lr
    iMant = 0
    MDCTLines1=np.zeros(codingParams.nMDCTLines)
    MDCTLines2=np.zeros(codingParams.nMDCTLines)

    for iBand in range(codingParams.sfBands.nBands):
        nLines = codingParams.sfBands.nLines[iBand]

        if iBand<codingParams.hBand and lr[iBand]==0:
            MDCTLines1[iMant:(iMant+nLines)] = (PreMDCTLine[0][iMant:(iMant+nLines)]+ PreMDCTLine[1][iMant:(iMant+ nLines)])
            MDCTLines2[iMant:(iMant+nLines)] = (PreMDCTLine[0][iMant:(iMant+nLines)]- PreMDCTLine[1][iMant:(iMant+ nLines)])
        else:
            MDCTLines1[iMant:(iMant+nLines)]= (PreMDCTLine[0][iMant:(iMant+nLines)])
            MDCTLines2[iMant:(iMant+nLines)]= (PreMDCTLine[1][iMant:(iMant+nLines)])

        iMant += nLines

    # IMDCT and window the data for this channel
    MDCTLines1 /= rescaleLevel  # put overall gain back to original level
    MDCTLines2 /= rescaleLevel  # put overall gain back to original level
    data1 = wn.SineWindow( IMDCT(MDCTLines1, halfN, halfN) )  # takes in halfN MDCT coeffs
    data2 = wn.SineWindow( IMDCT(MDCTLines2, halfN, halfN) )  # takes in halfN MDCT coeffs
    data = (data1,data2)

    # end loop over channels, return reconstituted time samples (pre-overlap-and-add)
    return data

def EncodeStereo(data,codingParams,remainder):

    # data to be appended to each channel
    ch1 = np.zeros(codingParams.nMDCTLines)
    ch2 = np.zeros(codingParams.nMDCTLines)
    nsamp=codingParams.nMDCTLines
    # number of lines in each band
    nLines = psy.AssignMDCTLinesFromFreqLimits(nsamp,codingParams.sampleRate)
    scaleFactorBands = psy.ScaleFactorBands(nLines)

    # MDCT of left and right channels
    lMDCT = MDCT(wn.SineWindow(data[0]),nsamp,nsamp)
    rMDCT = MDCT(wn.SineWindow(data[1]),nsamp,nsamp)

    #  generate masker SPL for frequencies specified
    # data: time domain samples
    # fft: windowed FFT function
    # freqs: frequencies lines desired for the mask values (Hz)
    # codingParams.sampleRate: sampling frequency (Hz)
    freqs = np.float_((np.arange(codingParams.nMDCTLines)+.5)*codingParams.sampleRate)/np.float(codingParams.nMDCTLines*2)  # MDCT frequency lines
    lMask = psy.Maskingcurve(data[0], freqs, codingParams.sampleRate)
    rMask = psy.Maskingcurve(data[1], freqs, codingParams.sampleRate)
    # decide whether critical band is L/R or M/S
    # compare the maximum value of the masking thresholds in each band up to hBand, inclusive
    lr = np.ones(codingParams.hBand) #self twenty five long, will later truncate
    for iBand in np.arange(codingParams.sfBands.nBands):
        dif = np.max(np.abs(lMask[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)]-rMask[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)]))
        if dif>=2 or iBand >= codingParams.hBand: #replace back later 031213
        # if the masking threshold difference is larger than 2, do l/r coding
            ch1[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)] = lMDCT[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)]
            ch2[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)] = rMDCT[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)]
        else:
            lr[iBand] = 0
            # if the masking threshold difference is smaller than 2, do m/s coding
            ch1[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)] = (lMDCT[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)]+rMDCT[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)])/2.
            ch2[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)] = (lMDCT[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)]-rMDCT[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)])/2.

    codingParams.lr = lr

    # compute THRm, THRs, and MLDm, MLD (if applicable)
    mMask = np.zeros(len(nLines))
    sMask = np.zeros(len(nLines))
    if np.sum(lr)-len(lr)<0: # (if any lr is 0)
        THRm =psy.Maskingcurve((lMDCT + rMDCT)/2., freqs, codingParams.sampleRate)
        THRs =psy.Maskingcurve((lMDCT-rMDCT)/2., freqs, codingParams.sampleRate)

        SSSm = psy.SpreadSignalEnergy((lMDCT + rMDCT)/2., THRm, freqs, codingParams.sampleRate)
        SSSs = psy.SpreadSignalEnergy((lMDCT - rMDCT)/2., THRs, freqs, codingParams.sampleRate)

        MLDF = psy.spreadsignalfactor(freqs)
        MLDm = SSSm+MLDF

        MLDs = SSSs+MLDF

        mMask = np.maximum(THRm, np.minimum(THRs,MLDs))
        sMask = np.maximum(THRs, np.minimum(THRm,MLDm))

    # compute final masking thresholds
    totalMask1 = np.zeros(len(lMDCT))
    totalMask2 = np.zeros(len(lMDCT))
    for iBand in np.arange(len(nLines)):
        if iBand >= codingParams.hBand or lr[iBand]==1:
            totalMask1[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)] = lMask[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)]
            totalMask2[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)] = rMask[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)]
        else:
            totalMask1[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)] = mMask[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)]
            totalMask2[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)] = sMask[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)]

    SPL1 = SPLfromMDCT(ch1)
    SPL2 = SPLfromMDCT(ch2)

    preSMR1 = SPL1 - totalMask1
    preSMR2 = SPL2 - totalMask2

    SMR1 = np.zeros(len(nLines))
    SMR2 = np.zeros(len(nLines))
    llim = scaleFactorBands.lowerLine
    ulim = scaleFactorBands.upperLine

    for ii in range(scaleFactorBands.nBands):
        idxlower=llim[ii]
        idxhigh=ulim[ii]
        SMR1[ii]=np.max(preSMR1[idxlower:idxhigh+1])
        SMR2[ii]=np.max(preSMR2[idxlower:idxhigh+1])

    # allocate bits (altogether, before dividing the bit allocation scheme into 2 channels)
    bitBudget = codingParams.targetBitsPerSample*nsamp-codingParams.hBand
    bitBudget -= codingParams.nScaleBits*(codingParams.sfBands.nBands+1)
    bitBudget -= codingParams.nMantSizeBits*codingParams.sfBands.nBands
    maxMantBits = 16

    bitAlloc,remainder = BitAlloc(bitBudget*2 + remainder, maxMantBits, 2*len(nLines), np.concatenate([nLines,nLines]), np.concatenate([SMR1,SMR2]))

    usedBits = np.sum(bitAlloc*np.concatenate([nLines,nLines]))
    bitAlloc = np.reshape(bitAlloc,[2,len(bitAlloc)/2])

    # compute scale factor
    scaleFactor = np.zeros((2,len(nLines)))
    for iBand in np.arange(len(nLines)):
        if iBand <= codingParams.hBand:
            scaleFactor[0][iBand] = (ScaleFactor(np.max(np.abs(ch1[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)])),codingParams.nScaleBits,bitAlloc[0][iBand]))
            scaleFactor[1][iBand] = (ScaleFactor(np.max(np.abs(ch2[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)])),codingParams.nScaleBits,bitAlloc[1][iBand]))
        else:
            scaleFactor[0][iBand] = (ScaleFactor(np.max(np.abs(ch1[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)])),codingParams.nScaleBits,bitAlloc[0][iBand]))
            scaleFactor[1][iBand] = (ScaleFactor(np.max(np.abs(ch2[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)])),codingParams.nScaleBits,bitAlloc[1][iBand]))

    # compute mantissa
    mantissa = np.zeros((2,len(lMDCT)))
    for iBand in np.arange(len(nLines)):
            mantissa[0][scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)] = vMantissa(ch1[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)],scaleFactor[0][iBand],codingParams.nScaleBits,bitAlloc[0][iBand])
            mantissa[1][scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)] = vMantissa(ch2[scaleFactorBands.lowerLine[iBand]:(scaleFactorBands.upperLine[iBand]+1)],scaleFactor[1][iBand],codingParams.nScaleBits,bitAlloc[1][iBand])
    overallScaleFactor = 0

    return (scaleFactor,bitAlloc,mantissa,overallScaleFactor,remainder,codingParams)

def dotProd(a,b): return np.sum(a*b)

def SPLfromMDCT(data):
        N = 2*len(data)
        intensity = np.maximum(4*np.square(data), 0.0000000000000001)
        return 96. + 10.* np.log10 (intensity)

def MDCTfromSPL(data): return np.sqrt((np.power(10.,data-96.)/10.)/4.)
