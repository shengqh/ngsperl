import csv

class GTExItem:
  def __init__(self, row):
    parts = row["snp"].split('_')
    self.Locus = parts[0] + "_" + parts[1]
    self.RefAllele = parts[2];
    self.AltAllele = parts[3];
    self.Gene = row["gene"].split('.')[0];
    self.Beta = float(row["beta"]);
    self.Key = self.Locus + ":" + self.Gene
    self.RefAltKey = self.Locus + ":" + self.RefAllele + ":" + self.AltAllele + ":" + self.Gene 
    self.AltRefKey = self.Locus + ":" + self.AltAllele + ":" + self.RefAllele + ":" + self.Gene 

def readGTExResult(fileName):
  result = list()
  with open(fileName, "r") as f:
    mycsv = csv.DictReader(f, delimiter="\t")
    for row in mycsv:
      item = GTExItem(row)
      result.append(item)
  return(result)
  
