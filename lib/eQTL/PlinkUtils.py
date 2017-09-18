class PlinkSNPItem:
  def __init__(self, row):
    parts = row.split('\t')
    self.Chromosome = parts[0]
    self.Name = parts[1]
    self.Position = parts[3];
    self.MinorAllele = parts[4]
    self.MajorAllele = parts[5]
    self.Locus = self.Chromosome + "_" + self.Position
    self.MinorKey = self.Locus + "_" + self.MinorAllele + "_" + self.MajorAllele
    self.MajorKey = self.Locus + "_" + self.MajorAllele + "_" + self.MinorAllele

def readPlinkSNP(bimFile):
  result = list()
  with open(bimFile, "r") as f:
    for row in f:
      item = PlinkSNPItem(row.rstrip())
      result.append(item)
  return(result)
  