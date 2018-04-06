import unittest
from SequenceUtils import calcIdenticalRate

class RolePrimerCriteria:
  def __init__(self, primerSequenceStart, primerSequence, primerSilimarRatio):
    self.PrimerSequence = primerSequence
    self.PrimerSequenceStart = primerSequenceStart
    self.PrimerSequenceEnd = self.PrimerSequenceStart + len(self.PrimerSequence)
    self.PrimerSimilarRatio = primerSilimarRatio
    
  def accept(self, sequence):
    curPrimerSequence = sequence[self.PrimerSequenceStart:self.PrimerSequenceStart + len(self.PrimerSequence)]
    similarRatio = calcIdenticalRate(self.PrimerSequence, curPrimerSequence)
    return (similarRatio >= self.PrimerSimilarRatio)
  
class TestRolePrimerCriteria(unittest.TestCase):
    def setUp(self):
        self.criteria = RolePrimerCriteria(6, "GCTGGACCGGTATGTAGCA", 0.9)
        self.sequence = "XXXXXXGCTGGACCGGTATGTAGCAGTACGTAGTCCTATTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAAGATTGCCATCGTTTGGGCAATATCAATAGGAGTTTCAGTTCCTATCCCTGTGATTGGACTGAGGGACGAAAGCAAAGTGTTCGTGAATAATACTACCTGCGTGCTCAATGACCCGAACTTCGTTCTCATCGGGTCCTTCGTGGCATTCTTCATCCCGTTGACAATTATGGTGATCACCTACTTCTTAACGATCTACGTCCTACGCCGTCAAGCTTTGAT"
      
    def testAcceptPrimerNoMismatch(self):
        self.assertTrue(self.criteria.accept(self.sequence))
      
    def testAcceptPrimerOneMismatch(self):
        seq = self.sequence[:7] + 'a' + self.sequence[8:]
        self.assertTrue(self.criteria.accept(seq))
        
    def testAcceptPrimerTwoMismatch(self):
        seq = self.sequence[:7] + 'aa' + self.sequence[9:]
        self.assertFalse(self.criteria.accept(seq))
