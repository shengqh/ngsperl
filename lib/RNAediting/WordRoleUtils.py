from difflib import SequenceMatcher
import unittest

class WordRole:
#   def __init__(self):
#     self.SampleName = ""
#     self.Barcode = ""
#     self.PrimerSequenceStart = ""
#     self.PrimerSequence = ""
#     self.PrimerSimilarRatio = 0.9
#     self.IdenticalStart = 0
#     self.IdenticalSequence = ""
#     self.SiteAllows = []

  def __init__(self, sampleName, barcode, primerSequenceStart, primerSequence, primerSilimarRatio, identicalStart, identicalSequence, siteAllows):
    self.SampleName = sampleName
    self.Barcode = barcode
    self.PrimerSequenceStart = primerSequenceStart
    self.PrimerSequence = primerSequence
    self.PrimerSimilarRatio = primerSilimarRatio
    self.IdenticalStart = identicalStart
    self.IdenticalSequence = identicalSequence
    self.SiteAllows = siteAllows
  
  def getWord(self, sequence):
    curBarcode = sequence[0:len(self.Barcode)]
    
    if self.Barcode != curBarcode :
      return("")
    
    curPrimerSequence = sequence[self.PrimerSequenceStart:self.PrimerSequenceStart + len(self.PrimerSequence)]
    similarRatio = SequenceMatcher(None, self.PrimerSequence, curPrimerSequence).ratio()
    if similarRatio < self.PrimerSimilarRatio:
      return("")
    
    curSequence = sequence[self.IdenticalStart:self.IdenticalStart + len(self.IdenticalSequence)]
    word = ""
    for idx in range(0, len(self.IdenticalSequence)):
      refBase = self.IdenticalSequence[idx]
      if refBase == 'R':
        if curSequence[idx] not in self.SiteAllows:
          return("")
        word = word + curSequence[idx]
      elif curSequence[idx] != refBase:
        return("")
    
    return(word)
      
        
class TestWordRole(unittest.TestCase):
    def testGetWordTurnee(self):
        role = WordRole("MalikG", "ATATGA", 6, "GCTGGACCGGTATGTAGCA", 0.9, 25, "RTRCGTRRTCCTRTTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAA", ['A', 'G'])
        sequence = "ATATGAGCTGGACCGGTATGTAGCAGTACGTAGTCCTATTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAAGATTGCCATCGTTTGGGCAATATCAATAGGAGTTTCAGTTCCTATCCCTGTGATTGGACTGAGGGACGAAAGCAAAGTGTTCGTGAATAATACTACCTGCGTGCTCAATGACCCGAACTTCGTTCTCATCGGGTCCTTCGTGGCATTCTTCATCCCGTTGACAATTATGGTGATCACCTACTTCTTAACGATCTACGTCCTACGCCGTCAAGCTTTGAT"
        curWord = role.getWord(sequence)
        self.assertEqual(curWord, "GAAGA")

    def testGetWordKayla(self):
        role = WordRole("Kayla1", "AGATAC", 7, "TTAACCCTCACTAAAGGGATTCTCA", 0.9, 43, "GTGATAAGGTCAATGRGGAGATGTATATA", ['A', 'G'])
        sequence = "AGATACATTAACCCTCACTAAAGGGATTCTCAGGATGTCCTTCGTGATAAGGTCAATGAGGAGATGTATATAGAAAGGTTATTTGATCAATGGTACAACAGCTCCATGAACATCATCTGCACGTGGCTGACCCTATAGTGAGTCGTATTAAGATCGGAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAGAAAAAGAAGAATGGAATGAAGACAATGAAAAGAAAAAAAAAAATAGAATAGTTAATAGTACATTACGGCGGAGGATAAGAGTAAGAGTAGAATGAATAGAG"
        curWord = role.getWord(sequence)
        self.assertEqual(curWord, "A")

    def testGetWordHussain(self):
        role = WordRole("Hussain1", "AACCAT", 7, "TTAACCCTCACTAAAGGGAGCTGAT", 0.9, 200, "RTRCGTRRTCCTRTTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAA", ['A', 'G'])
        sequence = "AACCATATTAACCCTCACTAAAGGGAGCTGATATGCTGGTGGGACTACTTGTCATGCCCCTGTCTCTGCTTGCAATTCTTTATGATTATGTCTGGCCTTTACCTAGATATTTGTGCCCCGTCTGGATTTCACTAGATGTGCTATTTTCAACTGCGTCCATCATGCACCTCTGCGCCATATCGCTGGACCGGTATGTAGCAGTGCGTAATCCTGTTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAAGATTGCCATCGTTTGGGCAATATCAATAGGAGTTTCAGTTCCTA"
        curWord = role.getWord(sequence)
        self.assertEqual(curWord, "GGAAG")
