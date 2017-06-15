import argparse
import sys
import logging
import os

def findSample(rnaseqFile, famFile, outputFile, logger):
  ids = []
  logger.info("reading %s ..." % (rnaseqFile))
  with open(rnaseqFile, 'r') as f:
    ids=f.readline().strip().split('\t')[2:]

  idmap = dict((v, 1) for v in ids)

  logger.info("reading %s ..." % (famFile))
  with open(famFile, 'r') as f:
    with open(outputFile, 'w') as wr:
      for line in f:
        parts = line.split(' ')
        if parts[0] in idmap:
          wr.write('%s %s\n' % (parts[0], parts[1]))

  logger.info("done.")
        
def main():
  parser = argparse.ArgumentParser(description="find common samples between RNASeq and SNP data",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  parser.add_argument('-r', '--rnaseq', action='store', nargs='?', help='RNASeq data file, only header will be used', required=True)
  parser.add_argument('-f', '--fam', action='store', nargs='?', help='Plink fam file', required=True)
  parser.add_argument('-o', '--output', action='store', nargs='?', help='Output file', required=True)
  
  args = parser.parse_args()
  
  rnaseqFile = args.rnaseq
  famFile=args.fam
  outputFile = args.output

  logger = logging.getLogger('findCommonSample')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  findSample(args.rnaseq, args.fam, args.output, logger)
  
if __name__ == "__main__":
    main()

