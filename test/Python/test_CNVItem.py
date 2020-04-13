import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../lib/CQS')))

from CNVItem import readCNVFile
import unittest

class BasicTestSuite(unittest.TestCase):
  def test_readCNVFile(self):
    cnvFile = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../data/cnv.txt'))
    cnvItems, samples = readCNVFile(cnvFile)
    self.assertEqual(9, len(cnvItems))
    ci1 = cnvItems[0]
    self.assertEqual("chr1", ci1.Chromosome)
    self.assertEqual("MMP23B", ci1.Gene)
    self.assertEqual(1633552, ci1.Start)
    self.assertEqual(1634150, ci1.End)
    self.assertEqual("chr1:1633697-1633962", ci1.getName())
    self.assertEqual({'P_273_03': 'DEL', 'P_273_13': 'DEL', 'P_273_15': 'DEL', 'P_273_39': 'DEL'}, ci1.SampleCNVMap)

    cnvItems, samples = readCNVFile(cnvFile, returnFullCNV=True)
    self.assertEqual(9, len(cnvItems))
    ci1 = cnvItems[0]
    self.assertEqual({'P_273_03': 'DEL,1,0,134', 'P_273_13': 'DEL,1,0,107', 'P_273_15': 'DEL,1,0,82', 'P_273_39': 'DEL,1,0,88'}, ci1.SampleCNVMap)

if __name__ == '__main__':
    unittest.main()