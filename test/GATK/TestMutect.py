import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../lib/GATK')))

from Mutect import MutectResult
import unittest
import logging

class BasicTestSuite(unittest.TestCase):
  def test_readFromFile_Mutect2_1(self):
    logger = logging.getLogger('test')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
    m2File = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../data/lesion_100.rmdup.indel.recal-filtered.vcf'))
    mr = MutectResult()
    mr.readFromFile(logger, "surface_100", m2File)
    self.assertEqual("lesion_100", mr.TumorSampleName)
    self.assertEqual("surface_100", mr.NormalSampleName)

  def test_readFromFile_Mutect2_2(self):
    logger = logging.getLogger('test')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
    m2File = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../data/lesion_100_2.rmdup.indel.recal-filtered.vcf'))
    mr = MutectResult()
    mr.readFromFile(logger, "crypts_100", m2File)
    self.assertEqual("lesion_100_2", mr.TumorSampleName)
    self.assertEqual("crypts_100", mr.NormalSampleName)

  def test_readFromFile_Mutect1_1(self):
    logger = logging.getLogger('test')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
    m2File = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../data/MAP.10711.somatic.vcf'))
    mr = MutectResult()
    mr.readFromFile(logger, "MAP.10711", m2File)
    self.assertEqual("MAP.10711_polyp", mr.TumorSampleName)
    self.assertEqual("MAP.10711_WBC", mr.NormalSampleName)

  def test_readFromFile_Mutect1_2(self):
    logger = logging.getLogger('test')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
    m2File = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../data/MAP.00696.somatic.vcf'))
    mr = MutectResult()
    mr.readFromFile(logger, "MAP.00696", m2File)
    self.assertEqual("MAP.00696_polyp", mr.TumorSampleName)
    self.assertEqual("MAP.00696_WBC", mr.NormalSampleName)

if __name__ == '__main__':
    unittest.main()