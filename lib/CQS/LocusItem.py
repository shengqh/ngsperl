class LocusItem(object):
  def __init__(self):
    self.Locus = ""
    self.Chromosome = ""
    self.Start = 0
    self.End = 0
    self._name = ""

  def setLocus(self, chromosome, start, end):
    self.Chromosome = chromosome
    self.Start = start
    self.End = end
    
  def getLocusString(self):
    return("%s:%d-%d" % (self.Chromosome, self.Start, self.End))
  
  def setLocusString(self, locus):
    #print(locus)
    self.Locus = locus
    parts = locus.split(":")
    #print(parts)
    self.Chromosome = parts[0]
    positions = parts[1].split("-")
    self.Start = int(positions[0])
    self.End = int(positions[1])

  def getName(self):
    if self._name == "":
      return self.getLocusString()
    else:
      return self._name

  def setName(self, name):
    self._name = name
    
  def getLocusFileString(self):
    return("%s_%d_%d" % (self.Chromosome, self.Start, self.End))
  
  def str(self):
    return self.getLocusString()

def readBedFile(fileName):
  """Read bed file to list of LocusItem
  
  Arguments:
      fileName {str} -- bed file name

  Returns:
      Array of LocusItem
  """
  result = []
  with open(fileName, "r") as fin:
    for line in fin:
      parts = line.rstrip().split('\t')
      chrom = parts[0]
      start = int(parts[1])
      end = int(parts[2])
            
      locus = LocusItem()
      locus.setLocus(chrom, start, end)

      if len(parts) > 4:
        locus.setName(parts[4])

      result.append(locus)

  return(result)

