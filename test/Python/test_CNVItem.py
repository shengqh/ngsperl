import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../lib/CQS')))

from CNVItem import readCNVFile
import unittest

class BasicTestSuite(unittest.TestCase):
  def test_getQueryMap(self):
    cnvFile = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../data/cnv.txt'))
    cnvItems = readCNVFile(cnvFile)
    self.assertEqual(3, len(cnvItems))

if __name__ == '__main__':
    unittest.main()