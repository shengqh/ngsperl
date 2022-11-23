import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../bin')))

from file_def import read_meta
import unittest

class BasicTestSuite(unittest.TestCase):
  def test_read_meta(self):
    meta_file = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../data/SraRunTable.txt'))
    meta_dic=read_meta(meta_file, 0, 25)
    print(meta_dic)

if __name__ == '__main__':
    unittest.main()