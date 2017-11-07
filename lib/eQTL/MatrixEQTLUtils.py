import csv
from PlinkUtils import readPlinkSNP

class MatrixEQTLItem:
  def __init__(self, row):
    self.SNP = row["SNP"]
    self.Gene = row["gene"]
    self.Beta = float(row["beta"])
    self.Tstat = float(row["t-stat"])
    self.Pvalue = float(row["p-value"])
    self.FDR = float(row["FDR"])
    self.Locus = '0_0'
    self.MajorAllele = ' '
    self.MinorAllele = ' '
    self.MajorKey = ""
    self.MinorKey = ""

def readMatrixEQTLResult(fileName):
  result = list()
  with open(fileName, "r") as f:
    mycsv = csv.DictReader(f, delimiter="\t")
    for row in mycsv:
      item = MatrixEQTLItem(row)
      result.append(item)
  return(result)
  
def fillMatrixEQTL(items, bimFile):
  plinkSnps = {p.Name:p for p in readPlinkSNP(bimFile)}
  for item in items:
    if not item.SNP in plinkSnps:
      raise Exception("Cannot find %s in bim file %s" % (item.SNP, bimFile))
    plinkItem = plinkSnps[item.SNP]
    item.MajorAllele = plinkItem.MajorAllele
    item.MinorAllele = plinkItem.MinorAllele
    item.Locus = plinkItem.Chromosome + "_" + str(plinkItem.Position)
    item.MajorKey = item.Locus + ":" + item.MajorAllele + ":" + item.MinorAllele + ":" + item.Gene 
    item.MinorKey = item.Locus + ":" + item.MinorAllele + ":" + item.MajorAllele + ":" + item.Gene 
