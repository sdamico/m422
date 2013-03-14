import numpy as np

SMR_THRESHOLD = -12

# Question 1.b)
def BitAllocUniform(bitBudget, maxMantBits, nBands, nLines, SMR):
    """
    Return a hard-coded vector that, in the case of the signal use in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are uniformely distributed for the mantissas.
    """
    return np.array([4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4])

def BitAllocConstSNR(bitBudget, maxMantBits, nBands, nLines, SMR):
    """
    Return a hard-coded vector that, in the case of the signal use in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are distributed for the mantissas to try and keep a constant
    quantization noise floor (assuming a noise floor 6 dB per bit below
    the peak SPL line in the scale factor band).
    """
    return (np.array([10, 11, 13, 16, 16, 16, 16, 16, 16, 11, 8, 7, 6, 5, 4, 4, 4, 13, 16,
                     4, 4, 16, 4, 4, 0])/2).astype(int)

def BitAllocConstMNR(bitBudget, maxMantBits, nBands, nLines, SMR):
    """
    Return a hard-coded vector that, in the case of the signal use in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are distributed for the mantissas to try and keep the quantization
    noise floor a constant distance below (or above, if bit starved) the
    masked threshold curve (assuming a quantization noise floor 6 dB per
    bit below the peak SPL line in the scale factor band).
    """
    return (np.array([11, 4, 5, 8, 10, 11, 10, 8, 7, 4, 4, 0, 0, 0, 0, 4, 4, 5, 4, 0, 4,
                     5, 4, 4, 8])/2).astype(int)

# Question 1.c)
def BitAlloc(bitBudget, maxMantBits, nBands, nLines, SMR):
    """
    Allocates bits to scale factor bands so as to flatten the NMR across the spectrum

       Arguments:
           bitBudget is total number of mantissa bits to allocate
           maxMantBits is max mantissa bits that can be allocated per line
           nBands is total number of scale factor bands
           nLines[nBands] is number of lines in each scale factor band
           SMR[nBands] is signal-to-mask ratio in each scale factor band

        Return:
            bits[nBands] is number of bits allocated to each scale factor band

        Logic:
           Maximizing SMR over blook gives optimization result that:
               R(i) = P/N + (1 bit/ 6 dB) * (SMR[i] - avgSMR)
           where P is the pool of bits for mantissas and N is number of bands
           This result needs to be adjusted if any R(i) goes below 2 (in which
           case we set R(i)=0) or if any R(i) goes above maxMantBits (in
           which case we set R(i)=maxMantBits).  (Note: 1 Mantissa bit is
           equivalent to 0 mantissa bits when you are using a midtread quantizer.)
           We will not bother to worry about slight variations in bit budget due
           rounding of the above equation to integer values of R(i).
    """
    bit_alloc = np.zeros(nBands, dtype=int)
    if len(SMR) == 0:
        return bit_alloc

    SMR[np.isneginf(SMR)] = -5000

    bands_filled = np.zeros(nBands, dtype=bool)
    #orig_bits = bitBudget
    for allocation_pass in range(2):
        while bitBudget > 0 and (np.min(bands_filled) == False):
            band_index = np.argmax(np.array([float(SMR[i]) if not bands_filled[i] else
                    np.finfo(float).min for i in range(len(SMR))]))
            if bitBudget >= nLines[band_index] and not bands_filled[band_index]:
                bitBudget -= nLines[band_index]
                SMR[band_index] -= 6.0
                bit_alloc[band_index]+=1
                if bit_alloc[band_index] == maxMantBits or SMR[band_index] < np.min(SMR)+SMR_THRESHOLD:
                    bands_filled[band_index] = True
            else:
                bands_filled[band_index] = True

        if allocation_pass == 0:
            bands_filled = (bit_alloc >= maxMantBits).astype(bool)

            bands_to_purge = bit_alloc == 1
            bands_zeroed = bit_alloc <= 1

            bands_filled[bands_zeroed] = True

            bit_alloc[bands_to_purge] = 0

            bitBudget += sum(nLines[bands_to_purge])


    #band = 0
    #totalbits = 0
    #for b in bit_alloc:
    #  totalbits += nLines[band]*b
    #  band += 1

    #print orig_bits - totalbits
    # disables dynamic bit rate
    #bitBudget = 0

    return bit_alloc, bitBudget

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":
    pass
