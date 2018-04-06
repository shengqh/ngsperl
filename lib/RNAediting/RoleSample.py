import unittest
from RolePrimerCriteria import RolePrimerCriteria 
from RoleSequenceCriteria import RoleSequenceCriteria

class RoleSample(object):
  def __init__(self, sampleName, barcode, primerCriteria, sequenceCriteria):
    self.SampleName = sampleName
    self.Barcode = barcode
    self.PrimerCriteria = primerCriteria
    self.SequenceCriteria = sequenceCriteria
    
  def acceptSample(self, sequence):
    curBarcode = sequence[0:len(self.Barcode)]
    if self.Barcode != curBarcode :
      return(False)

    return self.PrimerCriteria.accept(sequence)    
    
  def acceptSequence(self, sequence):
    return self.SequenceCriteria.accept(sequence)
  
  def getSequence(self, sequence):
    return self.SequenceCriteria.getSequence(sequence)   
  
class TestRoleSample(unittest.TestCase):
  def setUp(self):
    primerCriteria = RolePrimerCriteria(6, "GCTGGACCGGTATGTAGCA", 0.9)
    sequenceCriteria = RoleSequenceCriteria(25, "RTRCGTRRTCCTRTTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAA")
    self.role = RoleSample("MalikG", "ATATGA", primerCriteria, sequenceCriteria)
    
  def testAccept(self):
    sequence = "ATATGAGCTGGACCGGTATGTAGCAGTACGTAGTCCTATTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAAGATTGCCATCGTTTGGGCAATATCAATAGGAGTTTCAGTTCCTATCCCTGTGATTGGACTGAGGGACGAAAGCAAAGTGTTCGTGAATAATACTACCTGCGTGCTCAATGACCCGAACTTCGTTCTCATCGGGTCCTTCGTGGCATTCTTCATCCCGTTGACAATTATGGTGATCACCTACTTCTTAACGATCTACGTCCTACGCCGTCAAGCTTTGAT"
    self.assertTrue(self.role.acceptSample(sequence))
    sequence = "cTATGAGCTGGACCGGTATGTAGCAGTACGTAGTCCTATTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAAGATTGCCATCGTTTGGGCAATATCAATAGGAGTTTCAGTTCCTATCCCTGTGATTGGACTGAGGGACGAAAGCAAAGTGTTCGTGAATAATACTACCTGCGTGCTCAATGACCCGAACTTCGTTCTCATCGGGTCCTTCGTGGCATTCTTCATCCCGTTGACAATTATGGTGATCACCTACTTCTTAACGATCTACGTCCTACGCCGTCAAGCTTTGAT"
    self.assertFalse(self.role.acceptSample(sequence))
