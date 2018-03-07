import unittest
import collections
from SequenceUtils import calcIdenticalRate, getMismatchIndex

Word = collections.namedtuple('Word', ['word', 'discarded'])

class WordRole:
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
      return(Word("", False))
    
    curPrimerSequence = sequence[self.PrimerSequenceStart:self.PrimerSequenceStart + len(self.PrimerSequence)]
    similarRatio = calcIdenticalRate(self.PrimerSequence, curPrimerSequence)
    if similarRatio < self.PrimerSimilarRatio:
      return(Word("", False))
    
    curSequence = sequence[self.IdenticalStart:self.IdenticalStart + len(self.IdenticalSequence)]
    word = ""
    for idx in range(0, len(self.IdenticalSequence)):
      refBase = self.IdenticalSequence[idx]
      if refBase == 'R':
        if curSequence[idx] not in self.SiteAllows:
          return(Word("", True))
        word = word + curSequence[idx]
      elif curSequence[idx] != refBase:
        return(Word("", True))
    
    return((Word(word, False)))
        
class TestWordRole(unittest.TestCase):
    def testGetWordTurnee(self):
        role = WordRole("MalikG", "ATATGA", 6, "GCTGGACCGGTATGTAGCA", 0.9, 25, "RTRCGTRRTCCTRTTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAA", ['A', 'G'])
        sequence = "ATATGAGCTGGACCGGTATGTAGCAGTACGTAGTCCTATTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAAGATTGCCATCGTTTGGGCAATATCAATAGGAGTTTCAGTTCCTATCCCTGTGATTGGACTGAGGGACGAAAGCAAAGTGTTCGTGAATAATACTACCTGCGTGCTCAATGACCCGAACTTCGTTCTCATCGGGTCCTTCGTGGCATTCTTCATCCCGTTGACAATTATGGTGATCACCTACTTCTTAACGATCTACGTCCTACGCCGTCAAGCTTTGAT"
        curWord = role.getWord(sequence)
        self.assertEqual(curWord.word, "GAAGA")
        self.assertFalse(curWord.discarded)

    def testGetWordTurneeNotDiscard1(self):
        role = WordRole("MalikG", "ATATGA", 6, "GCTGGACCGGTATGTAGCA", 0.9, 25, "RTRCGTRRTCCTRTTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAA", ['A', 'G'])
        sequence = "TTATGAGCTGGACCGGTATGTAGCAGTACGTAGTCCTATTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAAGATTGCCATCGTTTGGGCAATATCAATAGGAGTTTCAGTTCCTATCCCTGTGATTGGACTGAGGGACGAAAGCAAAGTGTTCGTGAATAATACTACCTGCGTGCTCAATGACCCGAACTTCGTTCTCATCGGGTCCTTCGTGGCATTCTTCATCCCGTTGACAATTATGGTGATCACCTACTTCTTAACGATCTACGTCCTACGCCGTCAAGCTTTGAT"
        curWord = role.getWord(sequence)
        self.assertEqual(curWord.word, "")
        self.assertFalse(curWord.discarded)

    def testGetWordTurneeNotDiscard2(self):
        role = WordRole("MalikG", "ATATGA", 6, "GCTGGACCGGTATGTAGCA", 0.9, 25, "RTRCGTRRTCCTRTTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAA", ['A', 'G'])
        sequence = "ATATGAGCTaaACCGGTATGTAGCAGTACGTAGTCCTATTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAAGATTGCCATCGTTTGGGCAATATCAATAGGAGTTTCAGTTCCTATCCCTGTGATTGGACTGAGGGACGAAAGCAAAGTGTTCGTGAATAATACTACCTGCGTGCTCAATGACCCGAACTTCGTTCTCATCGGGTCCTTCGTGGCATTCTTCATCCCGTTGACAATTATGGTGATCACCTACTTCTTAACGATCTACGTCCTACGCCGTCAAGCTTTGAT"
        curWord = role.getWord(sequence)
        self.assertEqual(curWord.word, "")
        self.assertFalse(curWord.discarded)

    def testGetWordTurneeDiscard(self):
        role = WordRole("MalikG", "ATATGA", 6, "GCTGGACCGGTATGTAGCA", 0.9, 25, "RTRCGTRRTCCTRTTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAA", ['A', 'G'])
        sequence = "ATATGAGCTGGACCGGTATGTAGCAGTACGTAGaCCTATTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAAGATTGCCATCGTTTGGGCAATATCAATAGGAGTTTCAGTTCCTATCCCTGTGATTGGACTGAGGGACGAAAGCAAAGTGTTCGTGAATAATACTACCTGCGTGCTCAATGACCCGAACTTCGTTCTCATCGGGTCCTTCGTGGCATTCTTCATCCCGTTGACAATTATGGTGATCACCTACTTCTTAACGATCTACGTCCTACGCCGTCAAGCTTTGAT"
        curWord = role.getWord(sequence)
        self.assertEqual(curWord.word, "")
        self.assertTrue(curWord.discarded)

    def testGetWordKayla(self):
        role = WordRole("Kayla1", "AGATAC", 7, "TTAACCCTCACTAAAGGGATTCTCA", 0.9, 43, "GTGATAAGGTCAATGRGGAGATGTATATA", ['A', 'G'])
        sequence = "AGATACATTAACCCTCACTAAAGGGATTCTCAGGATGTCCTTCGTGATAAGGTCAATGAGGAGATGTATATAGAAAGGTTATTTGATCAATGGTACAACAGCTCCATGAACATCATCTGCACGTGGCTGACCCTATAGTGAGTCGTATTAAGATCGGAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAGAAAAAGAAGAATGGAATGAAGACAATGAAAAGAAAAAAAAAAATAGAATAGTTAATAGTACATTACGGCGGAGGATAAGAGTAAGAGTAGAATGAATAGAG"
        curWord = role.getWord(sequence)
        self.assertEqual(curWord.word, "A")
        self.assertFalse(curWord.discarded)

    def testGetWordHussain(self):
        role = WordRole("Hussain1", "AACCAT", 7, "TTAACCCTCACTAAAGGGAGCTGAT", 0.9, 200, "RTRCGTRRTCCTRTTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAA", ['A', 'G'])
        sequence = "AACCATATTAACCCTCACTAAAGGGAGCTGATATGCTGGTGGGACTACTTGTCATGCCCCTGTCTCTGCTTGCAATTCTTTATGATTATGTCTGGCCTTTACCTAGATATTTGTGCCCCGTCTGGATTTCACTAGATGTGCTATTTTCAACTGCGTCCATCATGCACCTCTGCGCCATATCGCTGGACCGGTATGTAGCAGTGCGTAATCCTGTTGAGCATAGCCGGTTCAATTCGCGGACTAAGGCCATCATGAAGATTGCCATCGTTTGGGCAATATCAATAGGAGTTTCAGTTCCTA"
        curWord = role.getWord(sequence)
        self.assertEqual(curWord.word, "GGAAG")
        self.assertFalse(curWord.discarded)
