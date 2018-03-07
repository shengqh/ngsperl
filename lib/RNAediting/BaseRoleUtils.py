from SequenceUtils import calcIdenticalRate, getMismatchIndex
import unittest
import collections
import math

BaseResult = collections.namedtuple('BaseResult', ['filled', 'discarded'])

class BaseRoleCriteria:
  def __init__(self, primerSequenceStart, primerSequence, primerSilimarRatio, identicalSequenceStart, identicalSequence, identicalSilimarRatio):
    self.PrimerSequence = primerSequence
    self.PrimerSequenceStart = primerSequenceStart
    self.PrimerSequenceEnd = self.PrimerSequenceStart + len(self.PrimerSequence)
    self.PrimerSimilarRatio = primerSilimarRatio
    self.IdenticalSequence = identicalSequence
    self.IdenticalSequenceStart = identicalSequenceStart
    self.IdenticalSequenceEnd = self.IdenticalSequenceStart + len(self.IdenticalSequence)
    self.IdenticalSilimarRatio = identicalSilimarRatio
  
class BaseRole:
  def __init__(self, sampleName, barcode, criteria):
    self.SampleName = sampleName
    self.Barcode = barcode
    self.Criteria = criteria
    self.BaseDict = {}
    self.TotalRead = 0
    self.DiscardRead = 0
  
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
      return(BaseResult(False, False))
    
    curPrimerSequence = sequence[self.Criteria.PrimerSequenceStart:self.Criteria.PrimerSequenceEnd]
    similarRatio = calcIdenticalRate(self.Criteria.PrimerSequence, curPrimerSequence)
    if similarRatio < self.Criteria.PrimerSimilarRatio:
      return(BaseResult(False, False))
    
    self.TotalRead = self.TotalRead + 1
    
    curIdenticalSequence = sequence[self.Criteria.IdenticalSequenceStart:self.Criteria.IdenticalSequenceEnd]
    mismatchIndex = getMismatchIndex(self.Criteria.IdenticalSequence, curIdenticalSequence)
    minLength = min(len(curIdenticalSequence), len(self.Criteria.IdenticalSequence))
    if 1 -(len(mismatchIndex) * 1.0 / minLength) < self.Criteria.IdenticalSilimarRatio:
      self.DiscardRead = self.DiscardRead + 1
      return(BaseResult(False, True))
    
    if len(mismatchIndex) == 0:  
      for idx in range(0, min(len(curIdenticalSequence), len(self.Criteria.IdenticalSequence))):
        curBase =curIdenticalSequence[idx]
        self.addBaseToDict(idx, curBase)
    else:
      for misIndex in mismatchIndex:
        curBase =curIdenticalSequence[misIndex]
        self.addBaseToDict(misIndex, curBase)
    
    return(BaseResult(True, False))
        
class TestBaseRole(unittest.TestCase):
    def testFillBaseDictKaylaNoMismatch(self):
        role = BaseRole("Kayla1", "AGATAC", BaseRoleCriteria(7, "TTAACCCTCACTAAAGGGATTCTCA", 0.9, 43, "GTGATAAGGTCAATGAGGAGATGTATATAGA", 0.9))
        sequence = "AGATACATTAACCCTCACTAAAGGGATTCTCAGGATGTCCTTCGTGATAAGGTCAATGAGGAGATGTATATAGAAAGGTTATTTGATCAATGGTACAACAGCTCCATGAACATCATCTGCACGTGGCTGACCCTATAGTGAGTCGTATTAAGATCGGAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAGAAAAAGAAGAATGGAATGAAGACAATGAAAAGAAAAAAAAAAATAGAATAGTTAATAGTACATTACGGCGGAGGATAAGAGTAAGAGTAGAATGAATAGAG"
        self.assertTrue(role.fillBaseDict(sequence).filled)
        self.assertEqual(len(role.Criteria.IdenticalSequence), len(role.BaseDict))
        for idx in range(0, len(role.Criteria.IdenticalSequence)):
          self.assertEqual(1, len(role.BaseDict[idx]))
          self.assertEqual(1, role.BaseDict[idx][role.Criteria.IdenticalSequence[idx]])
        #print(role.BaseDict)
        
    def testFillBaseDictKaylaThreeMismatch(self):
        role = BaseRole("Kayla1", "AGATAC", BaseRoleCriteria(7, "TTAACCCTCACTAAAGGGATTCTCA", 0.9, 43, "GTGATAAGGTCAATGAGGAGATGTATATAGA", 0.9))
        sequence = "AGATACATTAACCCTCACTAAAGGGATTCTCAGGATGTCCTTCactATAAGGTCAATGAGGAGATGTATATAGAAAGGTTATTTGATCAATGGTACAACAGCTCCATGAACATCATCTGCACGTGGCTGACCCTATAGTGAGTCGTATTAAGATCGGAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAGAAAAAGAAGAATGGAATGAAGACAATGAAAAGAAAAAAAAAAATAGAATAGTTAATAGTACATTACGGCGGAGGATAAGAGTAAGAGTAGAATGAATAGAG"
        res = role.fillBaseDict(sequence)
        self.assertTrue(res.filled)
        self.assertFalse(res.discarded)
        self.assertEqual(3, len(role.BaseDict))
        self.assertEqual(1, role.BaseDict[0]['a'])
        self.assertEqual(1, role.BaseDict[1]['c'])
        self.assertEqual(1, role.BaseDict[2]['t'])
        #print(role.BaseDict)
        
    def testFillBaseDictKaylaFourMismatch(self):
        role = BaseRole("Kayla1", "AGATAC", BaseRoleCriteria(7, "TTAACCCTCACTAAAGGGATTCTCA", 0.9, 43, "GTGATAAGGTCAATGAGGAGATGTATATAGA", 0.9))
        sequence = "AGATACATTAACCCTCACTAAAGGGATTCTCAGGATGTCCTTCacacATAAGGTCAATGAGGAGATGTATATAGAAAGGTTATTTGATCAATGGTACAACAGCTCCATGAACATCATCTGCACGTGGCTGACCCTATAGTGAGTCGTATTAAGATCGGAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAGAAAAAGAAGAATGGAATGAAGACAATGAAAAGAAAAAAAAAAATAGAATAGTTAATAGTACATTACGGCGGAGGATAAGAGTAAGAGTAGAATGAATAGAG"
        res = role.fillBaseDict(sequence)
        self.assertFalse(res.filled)
        self.assertTrue(res.discarded)
        self.assertEqual(0, len(role.BaseDict))
        #print(role.BaseDict)
