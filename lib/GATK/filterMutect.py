import argparse
import sys
import logging
import os
import errno
from asyncore import read

from Mutect import MultiMutectResult

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def filterMutect(logger, inputFile, outputFile, minNormalDepth, minTumorDepth, minMinorAlleleDepth): 
  check_file_exists(inputFile) 

  logger.info("Reading %s ..." % inputFile)

  mutect = MultiMutectResult()
  mutect.readFromFile(inputFile)

  with open(outputFile, "wt") as fout:
    for comment in mutect.Comments:
      fout.write("%s\n" % comment)
    
    for item in mutect.Data:
      passed = False
      for sample in item.Samples:
        if sample.NormalDepth >= minNormalDepth and sample.TumorDepth >= minTumorDepth and sample.minMinorAlleleDepth >= minMinorAlleleDepth:
          passed = True
          break

      if passed:
        fout.write(item.line)
    
def main():
  DEBUG=False
  NotDEBUG=not DEBUG

  parser = argparse.ArgumentParser(description="filter mutect result to keep tumor sample only.",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input vcf file', required=NotDEBUG)
  parser.add_argument('--min_normal_depth', action='store', type=int, default=8, help='Input minimum depth in normal sample')
  parser.add_argument('--min_tumor_depth', action='store', type=int, default=14, help='Input minimum depth in tumor sample')
  parser.add_argument('--min_minor_allele', action='store', type=int, default=3, help='Input minimum minor allele depth in tumor sample')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output vcf file", required=NotDEBUG)

  args = parser.parse_args()

  if DEBUG:
    args.input = "H:/shengquanhu/projects/20190610_Ciombior_ExomeSeq/combined.tumor.vcf"
    args.output = "H:/shengquanhu/projects/20190610_Ciombior_ExomeSeq/combined.filtred.vcf"

  logger = logging.getLogger('filterMutect')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  filterMutect(logger, args.input, args.output, args.min_normal_depth, args.min_tumor_depth, args.min_minor_allele)

  logger.info("done.")

if __name__ == "__main__":
    main()
