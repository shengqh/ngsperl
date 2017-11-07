from difflib import SequenceMatcher
import unittest

class BaseRole:
  def __init__(self, sampleName, barcode, primerSequenceStart, primerSequence, primerSilimarRatio, identicalStart, identicalSequence):
    self.SampleName = sampleName
    self.Barcode = barcode
    self.PrimerSequenceStart = primerSequenceStart
    self.PrimerSequence = primerSequence
    self.PrimerSimilarRatio = primerSilimarRatio
    self.IdenticalStart = identicalStart
    self.IdenticalSequence = identicalSequence
    self.BaseDict = {}
  
  def addBaseToDict(self, idx, curBase):
    if idx not in self.BaseDict:
      self.BaseDict[idx] = {}
    idxDic = self.BaseDict[idx]
      
    if curBase not in idxDic:
      idxDic[curBase] = 1
    else:
      idxDic[curBase] = idxDic[curBase] + 1
    
  def fillBaseDict(self, sequence):
    curBarcode = sequence[0:len(self.Barcode)]
    
    if self.Barcode != curBarcode :
      return(False)
    
    curPrimerSequence = sequence[self.PrimerSequenceStart:self.PrimerSequenceStart + len(self.PrimerSequence)]
    similarRatio = SequenceMatcher(None, self.PrimerSequence, curPrimerSequence).ratio()
    if similarRatio < self.PrimerSimilarRatio:
      return(False)
    
    curSequence = sequence[self.IdenticalStart:self.IdenticalStart + len(self.IdenticalSequence)]
    mismatched = 0
    mismatchedIndex = 0
    for idx in range(0, len(curSequence)):
      if self.IdenticalSequence[idx] != curSequence[idx]:
        mismatched = mismatched + 1
        mismatchedIndex = idx
    
    if mismatched > 1:
      return(False)
    
    if mismatched == 0:  
      for idx in range(0, len(self.IdenticalSequence)):
        curBase =curSequence[idx]
        self.addBaseToDict(idx, curBase)
    else:
      curBase =curSequence[mismatchedIndex]
      self.addBaseToDict(mismatchedIndex, curBase)
    
    return(True)
        
class TestBaseRole(unittest.TestCase):
    def testFillBaseDictKaylaNoMismatch(self):
        role = BaseRole("Kayla1", "AGATAC", 7, "TTAACCCTCACTAAAGGGATTCTCA", 0.9, 43, "GTGATAAGGTCAATGAGGAGATGTATATAGA")
        sequence = "AGATACATTAACCCTCACTAAAGGGATTCTCAGGATGTCCTTCGTGATAAGGTCAATGAGGAGATGTATATAGAAAGGTTATTTGATCAATGGTACAACAGCTCCATGAACATCATCTGCACGTGGCTGACCCTATAGTGAGTCGTATTAAGATCGGAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAGAAAAAGAAGAATGGAATGAAGACAATGAAAAGAAAAAAAAAAATAGAATAGTTAATAGTACATTACGGCGGAGGATAAGAGTAAGAGTAGAATGAATAGAG"
        self.assertTrue(role.fillBaseDict(sequence))
        self.assertEqual(len(role.IdenticalSequence), len(role.BaseDict))
        for idx in range(0, len(role.IdenticalSequence)):
          self.assertEqual(1, len(role.BaseDict[idx]))
          self.assertEqual(1, role.BaseDict[idx][role.IdenticalSequence[idx]])
        #print(role.BaseDict)
        
    def testFillBaseDictKaylaOneMismatch(self):
        role = BaseRole("Kayla1", "AGATAC", 7, "TTAACCCTCACTAAAGGGATTCTCA", 0.9, 43, "GTGATAAGGTCAATGAGGAGATGTATATAGA")
        sequence = "AGATACATTAACCCTCACTAAAGGGATTCTCAGGATGTCCTTCaTGATAAGGTCAATGAGGAGATGTATATAGAAAGGTTATTTGATCAATGGTACAACAGCTCCATGAACATCATCTGCACGTGGCTGACCCTATAGTGAGTCGTATTAAGATCGGAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAGAAAAAGAAGAATGGAATGAAGACAATGAAAAGAAAAAAAAAAATAGAATAGTTAATAGTACATTACGGCGGAGGATAAGAGTAAGAGTAGAATGAATAGAG"
        self.assertTrue(role.fillBaseDict(sequence))
        self.assertEqual(1, len(role.BaseDict))
        self.assertEqual(1, role.BaseDict[0]['a'])
        #print(role.BaseDict)
        
    def testFillBaseDictKaylaTwoMismatch(self):
        role = BaseRole("Kayla1", "AGATAC", 7, "TTAACCCTCACTAAAGGGATTCTCA", 0.9, 43, "GTGATAAGGTCAATGAGGAGATGTATATAGA")
        sequence = "AGATACATTAACCCTCACTAAAGGGATTCTCAGGATGTCCTTCaTaATAAGGTCAATGAGGAGATGTATATAGAAAGGTTATTTGATCAATGGTACAACAGCTCCATGAACATCATCTGCACGTGGCTGACCCTATAGTGAGTCGTATTAAGATCGGAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAGAAAAAGAAGAATGGAATGAAGACAATGAAAAGAAAAAAAAAAATAGAATAGTTAATAGTACATTACGGCGGAGGATAAGAGTAAGAGTAGAATGAATAGAG"
        self.assertFalse(role.fillBaseDict(sequence))
        self.assertEqual(0, len(role.BaseDict))
        #print(role.BaseDict)
