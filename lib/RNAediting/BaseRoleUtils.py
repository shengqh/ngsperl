from difflib import SequenceMatcher
import unittest

class BaseRole:
  def __init__(self, sampleName, barcode, primerSequenceStart, primerSequence, primerSilimarRatio, identicalStart, identicalSequence, identicalSimilarRatio, siteAllows):
    self.SampleName = sampleName
    self.Barcode = barcode
    self.PrimerSequenceStart = primerSequenceStart
    self.PrimerSequence = primerSequence
    self.PrimerSimilarRatio = primerSilimarRatio
    self.IdenticalStart = identicalStart
    self.IdenticalSequence = identicalSequence
    self.IdenticalSimilarRatio = primerSilimarRatio
    self.SiteAllows = siteAllows
    self.BaseDict = {}
  
  def fillBaseDict(self, sequence):
    curBarcode = sequence[0:len(self.Barcode)]
    
    if self.Barcode != curBarcode :
      return(False)
    
    curPrimerSequence = sequence[self.PrimerSequenceStart:self.PrimerSequenceStart + len(self.PrimerSequence)]
    similarRatio = SequenceMatcher(None, self.PrimerSequence, curPrimerSequence).ratio()
    if similarRatio < self.PrimerSimilarRatio:
      return(False)
    
    curSequence = sequence[self.IdenticalStart:self.IdenticalStart + len(self.IdenticalSequence)]
    similarRatio = SequenceMatcher(None, self.IdenticalSequence, curSequence).ratio()
    if similarRatio < self.PrimerSimilarRatio:
      return(False)
    
    for idx in range(0, len(self.IdenticalSequence)):
      refBase = self.IdenticalSequence[idx]
      curBase =curSequence[idx]
      if refBase == 'R':
        if curBase not in self.SiteAllows:
          return(False)
        
    for idx in range(0, len(self.IdenticalSequence)):
      refBase = self.IdenticalSequence[idx]
      curBase =curSequence[idx]

      if idx not in self.BaseDict:
        self.BaseDict[idx] = {}
      idxDic = self.BaseDict[idx]
      
      if curBase not in idxDic:
        idxDic[curBase] = 1
      else:
        idxDic[curBase] = idxDic[curBase] + 1
    
    return(True)
        
class TestBaseRole(unittest.TestCase):
    def testFillBaseDictKayla(self):
        role = BaseRole("Kayla1", "AGATAC", 7, "TTAACCCTCACTAAAGGGATTCTCA", 0.9, 43, "RTGATAAGGTCAATGAGGAGATGTATATAGA", 0.8, ['A', 'G'])
        baseDic = {}
        sequence = "AGATACATTAACCCTCACTAAAGGGATTCTCAGGATGTCCTTCGTGATAAGGTCAATGAGGAGATGTATATAGAAAGGTTATTTGATCAATGGTACAACAGCTCCATGAACATCATCTGCACGTGGCTGACCCTATAGTGAGTCGTATTAAGATCGGAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAGAAAAAGAAGAATGGAATGAAGACAATGAAAAGAAAAAAAAAAATAGAATAGTTAATAGTACATTACGGCGGAGGATAAGAGTAAGAGTAGAATGAATAGAG"
        self.assertTrue(role.fillBaseDict(sequence))
        print(role.BaseDict)
