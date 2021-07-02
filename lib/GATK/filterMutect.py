import argparse
import sys
import logging
import os
import errno
import gzip
from asyncore import read

from Mutect import MultiMutectItem

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def filterMutect(logger, inputFile, outputFile, minNormalDepth, minTumorDepth, minMinorAlleleDepth): 
  check_file_exists(inputFile) 

  logger.info("Processing %s ..." % inputFile)

  lineCount = 0
  passedCount = 0
  failedCount = 0

  tmpFile = outputFile + ".tmp.vcf"
  with open(tmpFile, "wt") as fout:
    if inputFile.endswith(".gz"):
      fin = gzip.open(inputFile,'rt')
    else:
      fin = open(inputFile, "rt")

    with fin:
      for line in fin:
        if line.startswith("#"):
          fout.write(line)
          continue

        lineCount += 1
        if lineCount % 10000 == 0:
          logger.info("%d" % lineCount)

        item = MultiMutectItem(line)

        passed = False
        for sample in item.Samples:
          if sample.NormalDepth == None:
            if sample.TumorDepth >= minTumorDepth and sample.MinorAlleleDepth >= minMinorAlleleDepth:
              passed = True
              break
          else:
            if sample.NormalDepth >= minNormalDepth and sample.TumorDepth >= minTumorDepth and sample.MinorAlleleDepth >= minMinorAlleleDepth:
              passed = True
              break

        if passed:
          fout.write(line)
          passedCount += 1
        else:
          failedCount += 1
      #  logger.info("Discarded:" + item.Line)
      
  with open(outputFile + ".stat", "wt") as fout:
    fout.write("minNormalDepth\t%d\n" % minNormalDepth)
    fout.write("minTumorDepth\t%d\n" % minTumorDepth)
    fout.write("minMinorAlleleDepth\t%d\n" % minMinorAlleleDepth)
    fout.write("Passed\t%d\n" % passedCount)
    fout.write("Failed\t%d\n" % failedCount)

  if os.path.isfile(outputFile):
    os.remove(outputFile)
  os.rename(tmpFile, outputFile)

  logger.info("Passed=%d, failed=%d" % (passedCount, failedCount))
    
def main():
  DEBUG=False
  NotDEBUG=not DEBUG

  parser = argparse.ArgumentParser(description="filter mutect result to keep tumor sample only.",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input vcf file', required=NotDEBUG)
  parser.add_argument('--min_normal_depth', action='store', type=int, default=8, help='Input minimum depth in normal sample')
  parser.add_argument('--min_tumor_depth', action='store', type=int, default=10, help='Input minimum depth in tumor sample')
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
