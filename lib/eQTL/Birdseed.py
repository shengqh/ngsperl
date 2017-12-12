import csv
import unittest

class BirdseedItem:
  def __init__(self, row):
    self.SNP = row["Probe Set ID"]
    self.SignalA = float(row["Signal A"])
    self.SignalB = float(row["Signal B"])
    self.RatioBA = self.SignalB / self.SignalA 
    calls = row["Forward Strand Base Calls"]
    x = calls[0]
    y = calls[1]
    if self.RatioBA < 1:
      if x == y:
        self.BaseA = x
        self.BaseB =  " "
        self.Genotype = 0
      else:
        self.BaseA = x
        self.BaseB = y
        self.Genotype = 1
    else:
      if x == y:
        self.BaseA = " "
        self.BaseB = y
        self.Genotype = 2
      else:
        self.BaseA = x
        self.BaseB = y
        self.Genotype = 1


def readBirdseed(fileName):
  result = list()
  fp = open(fileName)
  try:
    mycsv = csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter="\t")
    for row in mycsv:
      item = BirdseedItem(row)
      result.append(item)
  finally:
    fp.close()  
  return(result)

class TestBirdseedItem(unittest.TestCase):
  def testGenotype0(self):
    row = {"Probe Set ID":"SNP_A-2131660",
           "Forward Strand Base Calls":"CC",
           "Signal A":"2432.944",
           "Signal B":"487.201"}
    bi = BirdseedItem(row)
    self.assertEqual(bi.SNP, "SNP_A-2131660")
    self.assertEqual(bi.BaseA, "C")
    self.assertEqual(bi.BaseB, " ")
    self.assertEqual(bi.SignalA, 2432.944)
    self.assertEqual(bi.SignalB, 487.201)
    self.assertEqual(bi.Genotype, 0)
    
  def testGenotype1(self):
    row = {"Probe Set ID":"SNP_A-2131660",
           "Forward Strand Base Calls":"CT",
           "Signal A":"2432.944",
           "Signal B":"2487.201"}
    bi = BirdseedItem(row)
    self.assertEqual(bi.SNP, "SNP_A-2131660")
    self.assertEqual(bi.BaseA, "C")
    self.assertEqual(bi.BaseB, "T")
    self.assertEqual(bi.SignalA, 2432.944)
    self.assertEqual(bi.SignalB, 2487.201)
    self.assertEqual(bi.Genotype, 1)
    
  def testGenotype2(self):
    row = {"Probe Set ID":"SNP_A-2131660",
           "Forward Strand Base Calls":"TT",
           "Signal A":"432.944",
           "Signal B":"2487.201"}
    bi = BirdseedItem(row)
    self.assertEqual(bi.SNP, "SNP_A-2131660")
    self.assertEqual(bi.BaseA, " ")
    self.assertEqual(bi.BaseB, "T")
    self.assertEqual(bi.SignalA, 432.944)
    self.assertEqual(bi.SignalB, 2487.201)
    self.assertEqual(bi.Genotype, 2)
        