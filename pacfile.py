"""
pacfile.py -- Defines a PACFile class to handle reading and writing audio
data to an audio file holding data compressed using an MDCT-based perceptual audio
coding algorithm.  The MDCT lines of each audio channel are grouped into bands,
each sharing a single scaleFactor and bit allocation that are used to block-
floating point quantize those lines.  This class is a subclass of AudioFile.

-----------------------------------------------------------------------
© 2009 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------

See the documentation of the AudioFile class for general use of the AudioFile
class.

Notes on reading and decoding PAC files:

    The OpenFileForReading() function returns a CodedParams object containing:

        nChannels = the number of audio channels
        sampleRate = the sample rate of the audio samples
        numSamples = the total number of samples in the file for each channel
        nMDCTLines = half the MDCT block size (block switching not supported)
        nSamplesPerBlock = MDCTLines (but a name that PCM files look for)
        nScaleBits = the number of bits storing scale factors
        nMantSizeBits = the number of bits storing mantissa bit allocations
        sfBands = a ScaleFactorBands object
        overlapAndAdd = decoded data from the prior block (initially all zeros)

    The returned ScaleFactorBands object, sfBands, contains an allocation of
    the MDCT lines into groups that share a single scale factor and mantissa bit
    allocation.  sfBands has the following attributes available:

        nBands = the total number of scale factor bands
        nLines[iBand] = the number of MDCT lines in scale factor band iBand
        lowerLine[iBand] = the first MDCT line in scale factor band iBand
        upperLine[iBand] = the last MDCT line in scale factor band iBand


Notes on encoding and writing PAC files:

    When writing to a PACFile the CodingParams object passed to OpenForWriting()
    should have the following attributes set:

        nChannels = the number of audio channels
        sampleRate = the sample rate of the audio samples
        numSamples = the total number of samples in the file for each channel
        nMDCTLines = half the MDCT block size (format does not support block switching)
        nSamplesPerBlock = MDCTLines (but a name that PCM files look for)
        nScaleBits = the number of bits storing scale factors
        nMantSizeBits = the number of bits storing mantissa bit allocations
        targetBitsPerSample = the target encoding bit rate in units of bits per sample

    The first three attributes (nChannels, sampleRate, and numSamples) are
    typically added by the original data source (e.g. a PCMFile object) but
    numSamples may need to be extended to account for the MDCT coding delay of
    nMDCTLines and any zero-padding done in the final data block

    OpenForWriting() will add the following attributes to be used during the encoding
    process carried out in WriteDataBlock():

        sfBands = a ScaleFactorBands object
        priorBlock = the prior block of audio data (initially all zeros)

    The passed ScaleFactorBands object, sfBands, contains an allocation of
    the MDCT lines into groups that share a single scale factor and mantissa bit
    allocation.  sfBands has the following attributes available:

        nBands = the total number of scale factor bands
        nLines[iBand] = the number of MDCT lines in scale factor band iBand
        lowerLine[iBand] = the first MDCT line in scale factor band iBand
        upperLine[iBand] = the last MDCT line in scale factor band iBand

Description of the PAC File Format:

    Header:

        tag                 4 byte file tag equal to "PAC "
        sampleRate          little-endian unsigned long ("<L" format in struct)
        nChannels           little-endian unsigned short("<H" format in struct)
        numSamples          little-endian unsigned long ("<L" format in struct)
        nMDCTLines          little-endian unsigned long ("<L" format in struct)
        nScaleBits          little-endian unsigned short("<H" format in struct)
        nMantSizeBits       little-endian unsigned short("<H" format in struct)
        nSFBands            little-endian unsigned long ("<L" format in struct)
        for iBand in range(nSFBands):
            nLines[iBand]   little-endian unsigned short("<H" format in struct)

    Each Data Block:  (reads data blocks until end of file hit)

        for iCh in range(nChannels):
            nBytes          little-endian unsigned long ("<L" format in struct)
            as bits packed into an array of nBytes bytes:
                overallScale[iCh]                       nScaleBits bits
                for iBand in range(nSFBands):
                    scaleFactor[iCh][iBand]             nScaleBits bits
                    bitAlloc[iCh][iBand]                nMantSizeBits bits
                    if bitAlloc[iCh][iBand]:
                        for m in nLines[iBand]:
                            mantissa[iCh][iBand][m]     bitAlloc[iCh][iBand]+1 bits
                <extra custom data bits as long as space is included in nBytes>

"""

from audiofile import * # base class
from bitpack import *  # class for packing data into an array of bytes where each item's number of bits is specified
import codec as codec   # module where the actual PAC coding functions reside(this module only specifies the PAC file format)
from psychoac import ScaleFactorBands, AssignMDCTLinesFromFreqLimits  # defines the grouping of MDCT lines into scale factor bands
import sys
import huffman
import os
import numpy as np  # to allow conversion of data blocks to numpy's array object
MAX16BITS = 32767

class PACFile(AudioFile):
    """
    Handlers for a perceptually coded audio file I am encoding/decoding
    """

    # a file tag to recognize PAC coded files
    tag='PAC '

    def __init__(self, filename, training=False):
        AudioFile.__init__(self, filename)
        self.training = training

    def ReadFileHeader(self):
        """
        Reads the PAC file header from a just-opened PAC file and uses it to set
        object attributes.  File pointer ends at start of data portion.
        """

        # check file header tag to make sure it is the right kind of file
        tag=self.fp.read(4)
        if tag!=self.tag: raise "Tried to read a non-PAC file into a PACFile object"
        # use struct.unpack() to load up all the header data
        (sampleRate, nChannels, numSamples, nMDCTLines, nScaleBits, nMantSizeBits) \
                 = unpack('<LHLLHH',self.fp.read(calcsize('<LHLLHH')))
        nBands = unpack('<L',self.fp.read(calcsize('<L')))[0]
        nLines=  unpack('<'+str(nBands)+'H',self.fp.read(calcsize('<'+str(nBands)+'H')))
        sfBands=ScaleFactorBands(nLines)
        # load up a CodingParams object with the header data
        myParams=CodingParams()
        myParams.sampleRate = sampleRate
        myParams.nChannels = nChannels
        myParams.numSamples = numSamples
        myParams.nMDCTLines = myParams.nSamplesPerBlock = nMDCTLines
        myParams.nScaleBits = nScaleBits
        myParams.nMantSizeBits = nMantSizeBits
        # add in scale factor band information
        myParams.sfBands =sfBands
        # start w/o all zeroes as data from prior block to overlap-and-add for output
        overlapAndAdd = []
        for iCh in range(nChannels): overlapAndAdd.append( np.zeros(nMDCTLines, dtype=np.float64) )
        myParams.overlapAndAdd=overlapAndAdd
        return myParams

    def ReadDataBlock(self, codingParams):
        """
        Reads a block of coded data from a PACFile object that has already
        executed OpenForReading() and returns those samples as reconstituted
        signed-fraction data
        """
        # loop over channels (whose coded data are stored separately) and read in each data block
        data=[]
        numBands=25########hack added
        bitAllocfin=np.empty((2,25)) #######hack added
        mantfin=np.empty((2,codingParams.nMDCTLines))
        scalefactorfin=np.empty((2,25))
        codingParams.hBand = 16
        for iCh in range(codingParams.nChannels):
            data.append(np.array([],dtype=np.float64))  # add location for this channel's data
            # read in string containing the number of bytes of data for this channel (but check if at end of file!)
            s=self.fp.read(calcsize("<L"))  # will be empty if at end of file
            if not s:
                # hit last block, see if final overlap and add needs returning, else return nothing
                if codingParams.overlapAndAdd:
                    overlapAndAdd=codingParams.overlapAndAdd
                    codingParams.overlapAndAdd=0  # setting it to zero so next pass will just return
                    return overlapAndAdd
                else:
                    return
            # not at end of file, get nBytes from the string we just read
            nBytes = unpack("<L",s)[0] # read it as a little-endian unsigned long
            # read the nBytes of data into a PackedBits object to unpack
            pb = PackedBits()
            pb.SetPackedData( self.fp.read(nBytes) ) # PackedBits function SetPackedData() converts strings to internally-held array of bytes
            if pb.nBytes < nBytes:  raise "Only read a partial block of coded PACFile data"
            hdec = huffman.Huffman()
            lr=np.zeros(codingParams.sfBands.nBands)

            for iBand in range(codingParams.hBand): #Will need to change 25 to hBand
                lr[iBand]=pb.ReadBits(1)

            codingParams.lr=lr #store the information into codingParams
            ##############################this will work if hBand is set to maximum
            ####may need to add hband information

            hBand=pb.ReadBits(5)
            codingParams.hBand=hBand

            overallScaleFactor, scaleFactor, bitAlloc, mantissa = hdec.decodeMantissas(pb, codingParams)

            bitAllocfin[iCh]=bitAlloc #####hack added
            scalefactorfin[iCh]=scaleFactor
            mantfin[iCh]=mantissa

            # CUSTOM DATA:
            # < now can unpack any custom data passed in the nBytes of data >

            ##########################################adding extra info such as our lr data


            ######################################################

        # (DECODE HERE) decode the unpacked data for this channel, overlap-and-add first half, and append it to the data array (saving other half for next overlap-and-add)
        """plt.figure()
        plt.plot(mantfin[0],'g')
        plt.show()"""

        decodedData = self.Decode(scalefactorfin,bitAllocfin,mantfin, overallScaleFactor,codingParams)

        data[0]=np.concatenate( (data[0],np.add(codingParams.overlapAndAdd[0],decodedData[0][:codingParams.nMDCTLines]) ) )  # data[iCh] is overlap-and-added data
        data[1] = np.concatenate( (data[1],np.add(codingParams.overlapAndAdd[1],decodedData[1][:codingParams.nMDCTLines]) ) )  # data[iCh] is overlap-and-added data

        codingParams.overlapAndAdd[0] = decodedData[0][codingParams.nMDCTLines:]  # save other half for next pass
        codingParams.overlapAndAdd[1] = decodedData[1][codingParams.nMDCTLines:]  # save other half for next pass

        # end loop over channels, return signed-fraction samples for this block
        return data

    def WriteFileHeader(self,codingParams):
        """
        Writes the PAC file header for a just-opened PAC file and uses codingParams
        attributes for the header data.  File pointer ends at start of data portion.
        """
        self.remainder = 0
        self.huffman_training = huffman.HuffmanTrainer()
        # write a header tag
        self.fp.write(self.tag)
        # make sure that the number of samples in the file is a multiple of the
        # number of MDCT half-blocksize, otherwise zero pad as needed
        if not codingParams.numSamples%codingParams.nMDCTLines:
            codingParams.numSamples += (codingParams.nMDCTLines
                        - codingParams.numSamples%codingParams.nMDCTLines) # zero padding for partial final PCM block
        # also add in the delay block for the second pass w/ the last half-block
        codingParams.numSamples+= codingParams.nMDCTLines  # due to the delay in processing the first samples on both sides of the MDCT block
        # write the coded file attributes
        self.fp.write(pack('<LHLLHH',
            codingParams.sampleRate, codingParams.nChannels,
            codingParams.numSamples, codingParams.nMDCTLines,
            codingParams.nScaleBits, codingParams.nMantSizeBits  ))
        # create a ScaleFactorBand object to be used by the encoding process and write its info to header
        sfBands=ScaleFactorBands( AssignMDCTLinesFromFreqLimits(codingParams.nMDCTLines,
                                                                codingParams.sampleRate)
                                )
        codingParams.sfBands=sfBands
        self.fp.write(pack('<L',sfBands.nBands))
        self.fp.write(pack('<'+str(sfBands.nBands)+'H',*(sfBands.nLines.tolist()) ))
        # start w/o all zeroes as prior block of unencoded data for other half of MDCT block
        priorBlock = []
        for iCh in range(codingParams.nChannels):
            priorBlock.append(np.zeros(codingParams.nMDCTLines,dtype=np.float64) )
        codingParams.priorBlock = priorBlock
        return

    def WriteDataBlock(self,data, codingParams):
        """
        Writes a block of signed-fraction data to a PACFile object that has
        already executed OpenForWriting()"""

        # combine this block of multi-channel data w/ the prior block's to prepare for MDCTs twice as long
        fullBlockData=[]
        for iCh in range(codingParams.nChannels):
            fullBlockData.append( np.concatenate( ( codingParams.priorBlock[iCh], data[iCh]) ) )
        codingParams.priorBlock = data  # current pass's data is next pass's prior block data
        codingParams.hBand=16  ###################set hBand to max
        # (ENCODE HERE) Encode the full block of multi=channel data
        (scaleFactor,bitAlloc,mantissa, overallScaleFactor, self.remainder,codingParams) = self.Encode(fullBlockData,codingParams, self.remainder)  # returns a tuple with all the block-specific info not in the file header
        bitBudget = codingParams.targetBitsPerSample*codingParams.nMDCTLines*2
        # for each channel, write the data to the output file
        for iCh in range(codingParams.nChannels):
            pb = PackedBits()
            pb.Size(1)
            for iBand in range(codingParams.hBand):
                pb.WriteBits(int(codingParams.lr[iBand]),1)
            pb.WriteBits(int(codingParams.hBand),5)


            # Huffman encode the data
            henc = huffman.Huffman()

            pb = henc.encodeMantissas(pb, mantissa[iCh].astype(int),
                    scaleFactor[iCh].astype(int), codingParams.sfBands, codingParams, bitAlloc[iCh], overallScaleFactor)
            bitBudget -= pb.nBytes * 8

            self.fp.write(pack("<L",int(pb.nBytes))) # stores size as a little-endian unsigned long

            # finally, write the data in this channel's PackedBits object to the output file
            self.fp.write(pb.GetPackedData())

        if bitBudget < 0 and self.training:
            self.huffman_training.trainBlock(mantissa[0].astype(int), codingParams.sfBands,
                   codingParams, bitAlloc[0])
            self.huffman_training.trainBlock(mantissa[1].astype(int), codingParams.sfBands,
                   codingParams, bitAlloc[1])

        self.remainder += max(0, bitBudget)

        # end loop over channels, done writing coded data for all channels
        return

    def Close(self,codingParams):
        """
        Flushes the last data block through the encoding process (if encoding)
        and closes the audio file
        """
        # determine if encoding or encoding and, if encoding, do last block
        if self.fp.mode == "wb":  # we are writing to the PACFile, must be encode
            # we are writing the coded file -- pass a block of zeros to move last data block to other side of MDCT block
            data = [ np.zeros(codingParams.nMDCTLines,dtype=np.float),
                     np.zeros(codingParams.nMDCTLines,dtype=np.float) ]
            self.WriteDataBlock(data, codingParams)
            if self.training:
                self.OutputHuffmanTable()

        self.fp.close()

    def Encode(self,data,codingParams,remainder):
        """
        Encodes multichannel audio data and returns a tuple containing
        the scale factors, mantissa bit allocations, quantized mantissas,
        and the overall scale factor for each channel.
        """
        #Passes encoding logic to the Encode function defined in the codec module
        return codec.EncodeStereo(data,codingParams,remainder)

    def Decode(self,scaleFactor,bitAlloc,mantissa, overallScaleFactor,codingParams):
        """
        Decodes a single audio channel of data based on the values of its scale factors,
        bit allocations, quantized mantissas, and overall scale factor.
        """
        #Passes decoding logic to the Decode function defined in the codec module
        return codec.DecodeStereo(scaleFactor,bitAlloc,mantissa,0,codingParams)

    def GetNumberedFilename(self, prefix):
        filename_index = 0
        filename_attempt = prefix+'-'+str(filename_index)+'.csv'
        while filename_attempt in os.listdir('.'):
            filename_index += 1
            filename_attempt = prefix+'-'+str(filename_index)+'.csv'
        return filename_attempt

    def OutputHuffmanTable(self):
        symbol_dict = self.huffman_training.symbol_counts_bigvalue
        # TODO: support tables with != 8 bit symbols
        for i in range(256):
          if not i in symbol_dict:
            symbol_dict[i] = 0
        symbol_counts = np.array(symbol_dict.items())
        sorted_symbols = symbol_counts[symbol_counts[:,1].argsort()]
        sorted_symbols = sorted_symbols[::-1]
        if self.huffman_training.total_symbols_bigvalue > 0:
            symbol_prob = [[float(y)/self.huffman_training.total_symbols_bigvalue, x] for x, y in sorted_symbols]
            tree = huffman.buildTree(symbol_prob)
            huffman.saveTree(tree, self.GetNumberedFilename('bigvalue'))

        symbol_dict = self.huffman_training.symbol_counts_count1
        for i in range(16):
          if not i in symbol_dict:
            symbol_dict[i] = 0
        symbol_counts = np.array(symbol_dict.items())
        if len(symbol_counts) > 0:
          sorted_symbols = symbol_counts[symbol_counts[:,1].argsort()]
          sorted_symbols = sorted_symbols[::-1]
          if self.huffman_training.total_symbols_count1 > 0:
              symbol_prob = [[float(y)/self.huffman_training.total_symbols_count1, x] for x, y in sorted_symbols]
              tree = huffman.buildTree(symbol_prob)
              huffman.saveTree(tree, self.GetNumberedFilename('count1'))



#-----------------------------------------------------------------------------

# Testing the full PAC coder (needs a file called "input.wav" in the code directory)
if __name__=="__main__":

    print "\nTesting the PAC coder (input.wav -> coded.pac -> output.wav):"
    from pcmfile import * # to get access to WAV file handling

    for Direction in ("Encode","Decode"):

        # create the audio file objects
        if Direction == "Encode":
            print "\n\tEncoding input PCM file...",
            inFile= PCMFile("gsp.wav")
            outFile = PACFile("coded.pac")
        else:
            print "\n\tDecoding coded PAC file...",
            inFile = PACFile("coded.pac")
            outFile= PCMFile("gsp_decoded.wav")

        # open input file
        codingParams=inFile.OpenForReading()  # (includes reading header)

        # pass parameters to the output file
        if Direction == "Encode":
            # set additional parameters that are needed for PAC file (beyond those set by the PCM file on open)
            codingParams.nMDCTLines = 1024
            codingParams.nScaleBits = 3
            codingParams.nMantSizeBits = 4
            codingParams.targetBitsPerSample = int(np.floor(2.5*128000./float(codingParams.sampleRate)))
            # tell the PCM file how large the block size is
            codingParams.nSamplesPerBlock = codingParams.nMDCTLines
        else: # "Decode"
            # set PCM parameters (the rest is same as set by PAC file on open)
            codingParams.bitsPerSample = 16
        # only difference is in setting up the output file parameters

        # open the output file
        outFile.OpenForWriting(codingParams) # (includes writing header)

        # Read the input file and pass its data to the output file to be written
        #if Direction=="Decode"

        while True:
            data=inFile.ReadDataBlock(codingParams)
            if not data: break  # we hit the end of the input file
            outFile.WriteDataBlock(data,codingParams)
            print ".",  # just to signal how far we've gotten to user
            sys.stdout.flush()
        # end loop over reading/writing the blocks

            # close the files
        inFile.Close(codingParams)
        outFile.Close(codingParams)
        # end of loop over Encode/Decode

    print "\nDone with Encode/Decode test\n"
