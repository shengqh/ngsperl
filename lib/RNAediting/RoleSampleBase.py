import unittest
import collections
import math
from SequenceUtils import getMismatchIndex
from RoleSample import RoleSample
from RolePrimerCriteria import RolePrimerCriteria
from RoleSequenceCriteria import RoleSequenceCriteria

BaseResult = collections.namedtuple('BaseResult', ['filled', 'discarded'])

class RoleSampleBase(RoleSample):
  def __init__(self, sampleName, barcode, perimerCriteria, sequenceCriteria, identicalSilimarRatio):
    super( RoleSampleBase, self ).__init__(sampleName, barcode, perimerCriteria, sequenceCriteria)
    self.BaseDict = {}
    self.TotalRead = 0
    self.DiscardRead = 0
    self.IdenticalSilimarRatio = identicalSilimarRatio
  
  def addBaseToDict(self, idx, curBase):
    if idx not in self.BaseDict:
      self.BaseDict[idx] = {}
    idxDic = self.BaseDict[idx]
      
    if curBase not in idxDic:
      idxDic[curBase] = 1
    else:
      idxDic[curBase] = idxDic[curBase] + 1
    
  def fillBaseDict(self, sequence):
    if not self.acceptSample(sequence):
      return(BaseResult(False, False))
    
    self.TotalRead = self.TotalRead + 1
    
    curIdenticalSequence = self.getSequence(sequence)
    mismatchIndex = getMismatchIndex(self.SequenceCriteria.IdenticalSequence, curIdenticalSequence)
    minLength = min(len(curIdenticalSequence), len(self.SequenceCriteria.IdenticalSequence))
    if 1 -(len(mismatchIndex) * 1.0 / minLength) < self.IdenticalSilimarRatio:
      self.DiscardRead = self.DiscardRead + 1
      return(BaseResult(False, True))
    
    if len(mismatchIndex) == 0:  
      for idx in range(0, min(len(curIdenticalSequence), len(self.SequenceCriteria.IdenticalSequence))):
        curBase =curIdenticalSequence[idx]
        self.addBaseToDict(idx, curBase)
    else:
      for misIndex in mismatchIndex:
        curBase =curIdenticalSequence[misIndex]
        self.addBaseToDict(misIndex, curBase)
    
    return(BaseResult(True, False))
        
class TestRoleSampleBase(unittest.TestCase):
  def setUp(self):
    self.role = RoleSampleBase("Kayla1", "AGATAC", RolePrimerCriteria(7, "TTAACCCTCACTAAAGGGATTCTCA", 0.9), RoleSequenceCriteria(43, "GTGATAAGGTCAATGAGGAGATGTATATAGA"), 0.9)
    
  def testFillBaseDictKaylaNoMismatch(self):
      sequence = "AGATACATTAACCCTCACTAAAGGGATTCTCAGGATGTCCTTCGTGATAAGGTCAATGAGGAGATGTATATAGAAAGGTTATTTGATCAATGGTACAACAGCTCCATGAACATCATCTGCACGTGGCTGACCCTATAGTGAGTCGTATTAAGATCGGAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAGAAAAAGAAGAATGGAATGAAGACAATGAAAAGAAAAAAAAAAATAGAATAGTTAATAGTACATTACGGCGGAGGATAAGAGTAAGAGTAGAATGAATAGAG"
      self.assertTrue(self.role.fillBaseDict(sequence).filled)
      self.assertEqual(len(self.role.SequenceCriteria.IdenticalSequence), len(self.role.BaseDict))
      for idx in range(0, len(self.role.SequenceCriteria.IdenticalSequence)):
        self.assertEqual(1, len(self.role.BaseDict[idx]))
        self.assertEqual(1, self.role.BaseDict[idx][self.role.SequenceCriteria.IdenticalSequence[idx]])
      
  def testFillBaseDictKaylaThreeMismatch(self):
      sequence = "AGATACATTAACCCTCACTAAAGGGATTCTCAGGATGTCCTTCactATAAGGTCAATGAGGAGATGTATATAGAAAGGTTATTTGATCAATGGTACAACAGCTCCATGAACATCATCTGCACGTGGCTGACCCTATAGTGAGTCGTATTAAGATCGGAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAGAAAAAGAAGAATGGAATGAAGACAATGAAAAGAAAAAAAAAAATAGAATAGTTAATAGTACATTACGGCGGAGGATAAGAGTAAGAGTAGAATGAATAGAG"
      res = self.role.fillBaseDict(sequence)
      self.assertTrue(res.filled)
      self.assertFalse(res.discarded)
      self.assertEqual(3, len(self.role.BaseDict))
      self.assertEqual(1, self.role.BaseDict[0]['a'])
      self.assertEqual(1, self.role.BaseDict[1]['c'])
      self.assertEqual(1, self.role.BaseDict[2]['t'])
      
  def testFillBaseDictKaylaFourMismatch(self):
      sequence = "AGATACATTAACCCTCACTAAAGGGATTCTCAGGATGTCCTTCacacATAAGGTCAATGAGGAGATGTATATAGAAAGGTTATTTGATCAATGGTACAACAGCTCCATGAACATCATCTGCACGTGGCTGACCCTATAGTGAGTCGTATTAAGATCGGAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAGAAAAAGAAGAATGGAATGAAGACAATGAAAAGAAAAAAAAAAATAGAATAGTTAATAGTACATTACGGCGGAGGATAAGAGTAAGAGTAGAATGAATAGAG"
      res = self.role.fillBaseDict(sequence)
      self.assertFalse(res.filled)
      self.assertTrue(res.discarded)
      self.assertEqual(0, len(self.role.BaseDict))
