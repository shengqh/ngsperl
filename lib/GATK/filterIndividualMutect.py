import argparse
import sys
import logging
import os
import errno
import gzip

from Mutect import MutectResult

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def filterMutect(logger, fileName, filePath, outputFile, minNormalDepth, minTumorDepth, minAltAlleleDepth, minAltAlleleFrequency, filterFisherTest): 
  check_file_exists(filePath) 

  logger.info("Processing %s : %s ..." % (fileName, filePath))

  mutect = MutectResult()
  mutect.readFromFile(logger, fileName, filePath)  
  comments = mutect.Comments

  lineCount = 0
  passedCount = 0
  failedCount = 0
  fisherTestCount = 0

  tmpFile = outputFile + ".tmp.vcf"
  with open(tmpFile, "wt") as fout:
    for items in mutect.ChromosomeItemMap.values():
      for item in items:
        passed = True

        if item.TumorDepth < minTumorDepth or item.TumorAltAlleleDepth < minAltAlleleDepth or item.TumorAltAlleleFrequency < minAltAlleleFrequency:
          passed = False
        elif item.NormalDepth is not None:
          if item.NormalDepth < minNormalDepth:
            passed = False
          elif filterFisherTest:
            fpvalue = item.getFisherPValue()
            if fpvalue is not None and fpvalue > 0.05:
              fisherTestCount += 1
              passed = False

        if passed:
          fout.write(item.Line)
          passedCount += 1
        else:
          failedCount += 1
      
  with open(outputFile + ".stat", "wt") as fout:
    fout.write("minNormalDepth\t%d\n" % minNormalDepth)
    fout.write("minTumorDepth\t%d\n" % minTumorDepth)
    fout.write("minAltAlleleDepth\t%d\n" % minAltAlleleDepth)
    fout.write("minAltAlleleFrequency\t%f\n" % minAltAlleleFrequency)
    fout.write("Passed\t%d\n" % passedCount)
    fout.write("TotalFailed\t%d\n" % failedCount)
    fout.write("FisherTestFailed\t%d\n" % fisherTestCount)

  if os.path.isfile(outputFile):
    os.remove(outputFile)
  os.rename(tmpFile, outputFile)

  logger.info("Passed=%d, failed=%d, FisherTestFailed=%d" % (passedCount, failedCount, fisherTestCount))
    
def main():
  DEBUG=False
  NotDEBUG=not DEBUG

  parser = argparse.ArgumentParser(description="filter mutect result to keep tumor sample only.",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input vcf file', required=NotDEBUG)
  parser.add_argument('-n', '--sample_name', action='store', nargs='?', help='Input sample name', required=NotDEBUG)
  parser.add_argument('--min_normal_depth', action='store', type=int, default=8, help='Input minimum depth in normal sample')
  parser.add_argument('--min_tumor_depth', action='store', type=int, default=10, help='Input minimum depth in tumor sample')
  parser.add_argument('--min_alt_allele_depth', action='store', type=int, default=3, help='Input minimum alt allele depth in tumor sample')
  parser.add_argument('--min_alt_allele_frequency', action='store', type=float, default=0.05, help='Input minimum alt allele frequency in tumor sample')
  parser.add_argument('--filter_fisher_test', action='store', type=bool, default=False, help='Input whether to filter by Fisher\'s exact test')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output vcf file", required=NotDEBUG)

  args = parser.parse_args()

  if DEBUG:
    args.input = "H:/shengquanhu/projects/20190610_Ciombior_ExomeSeq/combined.tumor.vcf"
    args.output = "H:/shengquanhu/projects/20190610_Ciombior_ExomeSeq/combined.filtred.vcf"

  logger = logging.getLogger('filterMutect')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  filterMutect(logger, args.sample_name, args.input, args.output, args.min_normal_depth, args.min_tumor_depth, args.min_alt_allele_depth, args.min_alt_allele_frequency, args.filter_fisher_test)

  logger.info("done.")

if __name__ == "__main__":
    main()
