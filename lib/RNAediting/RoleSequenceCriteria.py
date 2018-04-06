import unittest
import collections
from SequenceUtils import calcIdenticalRate, getMismatchIndex

class RoleSequenceCriteria:
  def __init__(self, identicalSequenceStart, identicalSequence):
    self.IdenticalSequence = identicalSequence
    self.IdenticalSequenceStart = identicalSequenceStart
    self.IdenticalSequenceEnd = self.IdenticalSequenceStart + len(self.IdenticalSequence)
    
  def accept(self, sequence):
    return(True)
  
  def getSequence(self, sequence):
    return(sequence[self.IdenticalSequenceStart:self.IdenticalSequenceEnd])
        
class TestSequenceCriteria(unittest.TestCase):
    def setUp(self):
        self.criteria = RoleSequenceCriteria(25, "RTRCGTRRTCCTRTTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAA")
      
    def testGetSequence(self):
        sequence = "TTATGAGCTGGACCGGTATGTAGCAgTACGTAGTCCTATTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAAGATTGCCATCGTTTGGGCAATATCAATAGGAGTTTCAGTTCCTATCCCTGTGATTGGACTGAGGGACGAAAGCAAAGTGTTCGTGAATAATACTACCTGCGTGCTCAATGACCCGAACTTCGTTCTCATCGGGTCCTTCGTGGCATTCTTCATCCCGTTGACAATTATGGTGATCACCTACTTCTTAACGATCTACGTCCTACGCCGTCAAGCTTTGAT"
        self.assertEqual("gTACGTAGTCCTATTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAA", self.criteria.getSequence(sequence))
