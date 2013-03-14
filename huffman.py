from bitpack import *
import heapq
import csv
import copy
import os
BIGVALUE_MAX_BITS = 4
BIGVALUE_MAX = 15
COUNT1_MAX = 1
SYMBOL_MAX_LENGTH = 10
TABLE_ID_BITS = 3
REGION_BIGVALUE = 0
REGION_COUNT1 = 1
REGION_ZERO = 2
REGION_BITS = 5
LINBITS_SIZE_BITS = 4

# Inspired by this: http://en.literateprograms.org/index.php?title=Special:DownloadCode/Huffman_coding_(Python)&oldid=15672
def buildTree(probabilities_symbol_list):
  """Build a huffman tree given a list of symbols and probabilities.
  Inspired by this:
  http://en.literateprograms.org/index.php?title=Special:DownloadCode/
  Huffman_coding_(Python)&oldid=15672"""
  queue = list(probabilities_symbol_list)
  heapq.heapify(queue)
  while len(queue) > 1:
    right = heapq.heappop(queue)
    left = heapq.heappop(queue)
    parent = (right[0]+left[0], left, right)
    heapq.heappush(queue, parent)
  return queue[0]

def saveTree(tree, filename):
  with open(filename, 'w') as fp:
    printTree(tree, fp)

def printTree(tree, fp, prefix=''):
  """Output Huffman tree to file."""

  if len(tree) == 2:
    fp.write('%i,%s\n' %(tree[1], prefix))
  else:
    printTree(tree[1], fp, prefix + '0')
    printTree(tree[2], fp, prefix + '1')

def symbolBigValue(mantissas, bitAlloc):
  """Helper Function to look up the symbol for a big value mantissa"""

  mantissas_abs = mantissas & ~(1 << (bitAlloc-1))
  mantissas_capped = np.minimum(mantissas_abs, BIGVALUE_MAX)
  symbol = (mantissas_capped[1] << BIGVALUE_MAX_BITS) + mantissas_capped[0]
  return symbol

def symbolCount1(mantissas):
  """Helper Function to look up the symbol for a count1 mantissa"""

  mantissas_abs = np.minimum(mantissas.astype(int) & 1, 1)
  symbol = (mantissas_abs[3]<<3) + (mantissas_abs[2]<<2) + (mantissas_abs[1]<<1) + mantissas_abs[0]
  return symbol

def loadTable(filename):
  """Helper function to load a huffman table from a CSV."""

  table = {}
  with open(filename, 'rb') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
      table[int(row[0])] = (int(row[1],2), len(row[1]))
  return table


class HuffmanTrainer:
  """Class that performs huffman coding training."""

  def __init__(self):
    self.symbol_counts_bigvalue = {}
    self.symbol_counts_count1 = {}
    self.total_symbols_count1 = 0
    self.total_symbols_bigvalue = 0


  def trainBlock(self, mantissa, sfBands, codingParams, bitAlloc):
    iMant = 0
    max_bigvalue = 0
    count1_region = 0
    zero_region = 0
    for iBand in range(sfBands.nBands):
      if bitAlloc[iBand]:
        for j in range(sfBands.nLines[iBand]):
          mantissa_unsigned = mantissa[iMant+j] & ~(1 << (bitAlloc[iBand]-1))
          if mantissa_unsigned > 1:
            max_bigvalue = max(max_bigvalue, mantissa_unsigned)
            count1_region = iBand + 1
            zero_region = iBand + 1
          elif mantissa_unsigned > 0:
            zero_region = iBand + 1
      iMant += sfBands.nLines[iBand]
    linbits = 0
    if max_bigvalue > BIGVALUE_MAX:
      linbits = np.ceil(np.log2(max_bigvalue - BIGVALUE_MAX))

    iMant = 0
    for iBand in range(count1_region):
      if bitAlloc[iBand]:
        for j in range(0, sfBands.nLines[iBand], 2):
          mantissas = np.array([mantissa[iMant+j], 0])
          if j < (sfBands.nLines[iBand]-1):
            mantissas = np.array([mantissa[iMant+j], mantissa[iMant+j+1]])
          symbol = symbolBigValue(mantissas, bitAlloc[iBand])
          if symbol in self.symbol_counts_bigvalue:
            self.symbol_counts_bigvalue[symbol] += 1
          else:
            self.symbol_counts_bigvalue[symbol] = 1
          self.total_symbols_bigvalue += 1
      iMant = iMant + sfBands.nLines[iBand]

    for iBand in range(count1_region,zero_region):
      if bitAlloc[iBand]:
        for j in range(0, sfBands.nLines[iBand], 2):
          mantissas = np.array([mantissa[iMant+j], 0])
          num_mantissas = 4
          if j >= (sfBands.nLines[iBand]-3):
            num_mantissas = sfBands.nLines[iBand]-j

          mantissas = np.array(mantissa[iMant+j:iMant+j+num_mantissas])
          if num_mantissas < 4:
            mantissas = np.concatenate([mantissas, np.zeros(4-num_mantissas)])

          symbol = symbolCount1(mantissas)
          if symbol in self.symbol_counts_count1:
            self.symbol_counts_count1[symbol] += 1
          else:
            self.symbol_counts_count1[symbol] = 1
          self.total_symbols_count1 += 1
      iMant = iMant + sfBands.nLines[iBand]


class Huffman:
  def __init__(self):
    self.reverse_tables_bigvalue = []
    filenames = os.listdir('.')
    bigvalue_filenames = []
    count1_filenames = []
    for f in filenames:
      if "bigvalue" in f:
        bigvalue_filenames.append(f)
      elif "count1" in f:
        count1_filenames.append(f)

    bigvalue_filenames.sort()
    count1_filenames.sort()

    self.tables_bigvalue = [loadTable(f) for f in bigvalue_filenames]
    self.reverse_tables_count1 = []
    for table in self.tables_bigvalue:
      self.reverse_tables_bigvalue.append({(y, z): x for (x, (y, z)) in table.items()})
    self.tables_count1 = [loadTable(f) for f in count1_filenames]
    for table in self.tables_count1:
      self.reverse_tables_count1.append({(y, z): x for (x, (y, z)) in table.items()})

    self.bits_saved = 0

  def encodeMantissas(self, pb, mantissa, scaleFactor, sfBands, codingParams, bitAlloc, overallScaleFactor):
    iMant = 0
    max_bigvalue = 0
    count1_region = 0
    zero_region = 0
    mantissas_unsigned = []
    mantissas = []
    for iBand in range(sfBands.nBands):
      if bitAlloc[iBand]:
        for j in range(sfBands.nLines[iBand]):
          mantissa_unsigned = mantissa[iMant+j] & ~(1 << (bitAlloc[iBand]-1))
          mantissas.append(mantissa[iMant+j])
          mantissas_unsigned.append(mantissa_unsigned)
          if mantissa_unsigned > 1:
            max_bigvalue = max(max_bigvalue, mantissa_unsigned)
            count1_region = iBand + 1
            zero_region = iBand + 1
          elif mantissa_unsigned > 0:
            zero_region = iBand + 1
      iMant += sfBands.nLines[iBand]
    linbits = 0
    if max_bigvalue > BIGVALUE_MAX:
      linbits = int(np.ceil(np.log2(max_bigvalue - BIGVALUE_MAX)))
    bits_to_allocate = 0

    # Get bit allocation for big value region
    min_table = -1
    best_bitstream = None
    best_iMant = 0
    for cur_table in range(len(self.tables_bigvalue)):
      bitstream = copy.deepcopy(pb)
      bitstream.WriteBits(overallScaleFactor, codingParams.nScaleBits)
      bitstream.WriteBits(count1_region, int(REGION_BITS))
      bitstream.WriteBits(zero_region, int(REGION_BITS))
      bitstream.WriteBits(int(linbits), int(LINBITS_SIZE_BITS))
      bitstream.WriteBits(cur_table, TABLE_ID_BITS)
      iMant = 0
      for iBand in range(count1_region):
        ba = bitAlloc[iBand]
        if ba: ba-=1
        bitstream.WriteBits(ba, codingParams.nMantSizeBits)
        bitstream.WriteBits(scaleFactor[iBand], codingParams.nScaleBits)
        if bitAlloc[iBand]:
          for j in range(0, sfBands.nLines[iBand], 2):
            mantissas = np.array([mantissa[iMant+j], 0])
            if j < (sfBands.nLines[iBand]-1):
              mantissas = np.array([mantissa[iMant+j], mantissa[iMant+j+1]])
            self.encodeMantissasBigValue(mantissas, bitAlloc[iBand], bitstream, BIGVALUE_MAX, linbits, self.tables_bigvalue[cur_table])
        iMant = iMant + sfBands.nLines[iBand]
      if (min_table == -1) or bitstream.nBytes < best_bitstream.nBytes:
        min_table = cur_table
        best_bitstream = bitstream
        best_iMant = iMant

    # Get bit allocation for count1 region
    min_table_count1 = -1
    best_bitstream_count1 = None
    for cur_table in range(len(self.tables_count1)):
      bitstream = copy.deepcopy(best_bitstream)
      iMant = best_iMant
      bitstream.WriteBits(cur_table, TABLE_ID_BITS)
      for iBand in range(count1_region, zero_region):
        ba = bitAlloc[iBand]
        if ba: ba=1
        bitstream.WriteBits(ba, 1) #codingParams.nMantSizeBits)
        bitstream.WriteBits(scaleFactor[iBand], codingParams.nScaleBits)
        if bitAlloc[iBand]:
          for j in range(0, sfBands.nLines[iBand], 4):
            mantissas = np.array([mantissa[iMant+j], 0])
            num_mantissas = 4
            if j >= (sfBands.nLines[iBand]-3):
              num_mantissas = sfBands.nLines[iBand]-j

            mantissas = np.array(mantissa[iMant+j:iMant+j+num_mantissas])
            if num_mantissas < 4:
              mantissas = np.concatenate([mantissas, np.zeros(4-num_mantissas)])
            #print "encoding count1 mantissas", mantissas
            self.encodeMantissasCount1(mantissas, bitstream, self.tables_count1[cur_table])
        iMant = iMant + sfBands.nLines[iBand]
      if (min_table_count1 == -1) or bitstream.nBytes < best_bitstream.nBytes:
        min_table_count1 = cur_table
        best_bitstream_count1 = bitstream

    bitstream = best_bitstream_count1
    if best_bitstream_count1 == None:
      bitstream = best_bitstream

    return bitstream

  def decodeMantissas(self, bitstream, codingParams):
    overallScaleFactor = bitstream.ReadBits(codingParams.nScaleBits)  # overall scale factor
    count1_region = bitstream.ReadBits(int(REGION_BITS))
    zero_region = bitstream.ReadBits(int(REGION_BITS))

    linbits = bitstream.ReadBits(LINBITS_SIZE_BITS)
    scaleFactor=[]
    bitAlloc=[]
    mantissa=np.zeros(codingParams.nMDCTLines,np.int32)  # start w/ all mantissas zero
    cur_table = bitstream.ReadBits(TABLE_ID_BITS)


    for iBand in range(zero_region): # loop over each scale factor band to pack its data
      region = REGION_BIGVALUE
      if iBand >= count1_region:
        region = REGION_COUNT1
      if iBand >= zero_region:
        region = REGION_ZERO

      # Start of count1 region, read table ID
      if iBand == count1_region:
        cur_table = bitstream.ReadBits(TABLE_ID_BITS)


      if region == REGION_BIGVALUE:
        ba = bitstream.ReadBits(codingParams.nMantSizeBits)
        if ba: ba+=1  # no bit allocation of 1 so ba of 2 and up stored as one less
        bitAlloc.append(ba)  # bit allocation for this band
      elif region == REGION_COUNT1:
        ba = bitstream.ReadBits(1)
        if ba: ba+=1
        bitAlloc.append(ba)
      scaleFactor.append(bitstream.ReadBits(codingParams.nScaleBits))  # scale factor for this band

      if bitAlloc[iBand]:
        # if bits allocated, extract those mantissas and put in correct location in matnissa array
        m=np.empty(codingParams.sfBands.nLines[iBand],np.int32)
        if region == REGION_BIGVALUE:
          for j in range(0,codingParams.sfBands.nLines[iBand],2):
            #TODO: if there's an error and the symbol isn't found, we will
            # overflow our bitstream
            mantissas = self.decodeMantissasBigValue(bitstream, 32, linbits, bitAlloc[iBand], self.reverse_tables_bigvalue[cur_table])
            m[j]=mantissas[0]
            if (j+1) < codingParams.sfBands.nLines[iBand]:
              m[j+1]=mantissas[1]
        elif region == REGION_COUNT1:
          for j in range(0,codingParams.sfBands.nLines[iBand],4):
            mantissas = self.decodeMantissasCount1(bitstream, 32, self.reverse_tables_count1[cur_table])
            m[j]=mantissas[0]
            for n in range(1,4):
              if (j+n) < codingParams.sfBands.nLines[iBand]:
                m[j+n]=mantissas[n]
        mantissa[codingParams.sfBands.lowerLine[iBand]:(codingParams.sfBands.upperLine[iBand]+1)] = m

    # done unpacking data (end loop over scale factor bands)
    for iBand in range(zero_region, codingParams.sfBands.nBands):
      m=np.empty(codingParams.sfBands.nLines[iBand],np.int32)
      bitAlloc.append(0)
      scaleFactor.append(0)
      for j in range(codingParams.sfBands.nLines[iBand]):
        m[j]=0
      mantissa[codingParams.sfBands.lowerLine[iBand]:(codingParams.sfBands.upperLine[iBand]+1)] = m

    return overallScaleFactor, scaleFactor, bitAlloc, mantissa

  def encodeMantissasBigValue(self, mantissas, bitAlloc, bitstream, max_table_value, linbits, table):
    signs = mantissas >> (bitAlloc-1)
    mantissas_abs = (mantissas & ~(1 << (bitAlloc-1))).astype(int)
    mantissas_capped = (np.minimum(mantissas_abs, BIGVALUE_MAX)).astype(int)
    symbol = (mantissas_capped[1] << BIGVALUE_MAX_BITS) + mantissas_capped[0]
    linbits_values = mantissas_abs - mantissas_capped
    code, code_bits = table[symbol]
    bitstream.WriteBits(code, code_bits)
    for sign, mantissa_abs, mantissa_capped, linbits_value in zip(signs, mantissas_abs, mantissas_capped, linbits_values):
      if mantissa_capped == BIGVALUE_MAX:
        bitstream.WriteBits(linbits_value, int(linbits))
      if mantissa_abs > 0:
        bitstream.WriteBits(sign, 1)

  def encodeMantissasCount1(self, mantissas, bitstream, table):
    signs = mantissas.astype(int) >> int(1)
    mantissas_abs = np.minimum(mantissas.astype(int) & 1, 1)
    symbol = (mantissas_abs[3]<<3) + (mantissas_abs[2]<<2) + (mantissas_abs[1]<<1) + mantissas_abs[0]
    code, code_bits = table[symbol]
    bitstream.WriteBits(code, code_bits)
    for sign, mantissa_abs in zip(signs, mantissas_abs):
      if mantissa_abs > 0:
        bitstream.WriteBits(sign, 1)

  def decodeMantissasBigValue(self, bitstream, max_table_bits, linbits, bitAlloc, reverse_table):
    code = 0
    symbol = -1
    bits = 0
    for i in range(max_table_bits):
      bit = bitstream.ReadBits(1)
      code = (code<<1) | bit
      bits = i+1
      if (code, bits) in reverse_table:
        symbol = reverse_table[(code, bits)]
        break
    if symbol == -1:
      print "FAILURE FINDING SYMBOL"
      return None
    values = [symbol & ((1<<BIGVALUE_MAX_BITS)-1), symbol >> BIGVALUE_MAX_BITS]
    results = []
    for value in values:
      if value == 15 and linbits > 0:
        linbits_value = bitstream.ReadBits(linbits)
        value += linbits_value
      if value != 0:
        sign = bitstream.ReadBits(1)
        if sign:
          value = value + (1<<(bitAlloc-1))
      results.append(value)

    return results

  def decodeMantissasCount1(self, bitstream, max_table_bits, reverse_table):
    code = 0
    symbol = -1
    bits = 0
    for i in range(max_table_bits):
      bit = bitstream.ReadBits(1)
      code = (code<<1) | bit
      bits = i+1
      if (code, bits) in reverse_table:
        symbol = reverse_table[(code, bits)]
        break
    if symbol == -1:
      return None

    values = [symbol & 1, (symbol >> 1) & 1, (symbol >> 2) & 1, (symbol >> 3) & 1]

    results = []
    for value in values:
      if value != 0:
        sign = bitstream.ReadBits(1)
        if sign:
          value = value + (1<<1)
      results.append(value)

    return results
