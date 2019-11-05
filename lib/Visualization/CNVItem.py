from LocusItem import LocusItem

class CNVItem(LocusItem):
  def __init__(self, locus, name, sampleCNVMap):
    self.setLocusString(locus)
    self.Name = name
    self.SampleCNVMap = sampleCNVMap
    
  def overlap(self, locus):
    if self.Chromosome != locus.Chromosome:
      return(False)
    if self.End < locus.Start:
      return(False)
    if self.Start > locus.End:
      return(False)
    return(True)
  
  def contains(self, chromosome, position):
    if self.Chromosome != chromosome:
      return(False)
    if self.Start > position:
      return(False)
    if self.End < position:
      return(False)
    return(True)

