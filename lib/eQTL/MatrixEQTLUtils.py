import csv

class MatrixEQTLItem:
  def __init__(self, row):
    self.SNP = row["SNP"]
    self.Gene = row["gene"]
    self.Beta = float(row["beta"])
    self.Tstat = float(row["t-stat"])
    self.Pvalue = float(row["p-value"])
    self.FDR = float(row["FDR"])

def readMatrixEQTLResult(fileName):
  result = list()
  with open(fileName, "r") as f:
    mycsv = csv.DictReader(f, delimiter="\t")
    for row in mycsv:
      item = MatrixEQTLItem(row)
      result.append(item)
  return(result)
  