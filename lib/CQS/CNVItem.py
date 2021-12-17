from LocusItem import LocusItem

class CNVItem(LocusItem):
  def __init__(self):
    LocusItem.__init__(self)
    self.Gene = ""
    self.CytoBand = ""
    self.SampleCNVMap = {}

def readCNVFile(fileName, returnFullCNV=False):
  result = []
  with open(fileName, "r") as fin:
    headers = fin.readline().rstrip().split('\t')
    hasCytoBand = "cytoBand" == headers[5]
    startIndex = 6 if hasCytoBand else 5
    samples = headers[startIndex:]

    for line in fin:
      parts = line.rstrip().split('\t')
      chrom = parts[0]
      start = int(parts[1])
      end = int(parts[2])

      ci = CNVItem()
      ci.setLocus(chrom, start, end)
      ci.setName(parts[3])
      ci.Gene = parts[4]
      if hasCytoBand:
        ci.CytoBand = parts[5]

      sampleCNVMap = {}
      for sampleIndex in range(startIndex, len(parts)):
        cnv = parts[sampleIndex]
        if cnv != "":
          if returnFullCNV:
            sampleCNVMap[headers[sampleIndex]] = cnv
          else:
            sampleCNVMap[headers[sampleIndex]] = cnv.split(',', 1)[0]

      ci.SampleCNVMap = sampleCNVMap

      result.append(ci)

  return(result, samples)
