import sys
import os
import logging
import argparse
import string
import subprocess
import pysam

def count(logger, bam_file, output_file):
  with pysam.Samfile(bam_file) as samfile:
    with open(output_file, "w") as fout:
      fout.write("Chromosome\tCount\n")
      logger.info("start counting %s ..." % bam_file)
      curChromosomes = samfile.references
      for chromosome in curChromosomes:
        chromCount = samfile.count(chromosome)
        logger.info("%s : %d" % (chromosome, chromCount))
        fout.write("%s\t%d\n" % (chromosome, chromCount))
          
  logger.info("done.")

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Get read count in each chromosome.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help="Input bam file", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output count file", required=NOT_DEBUG)

  args = parser.parse_args()
  
  if(DEBUG):
    args.input = "/scratch/jbrown_lab/shengq2/projects/20200905_wgs_5162_hg38/bwa/result/CF001.sortedByCoord.bam"
    args.output = "/scratch/jbrown_lab/shengq2/projects/20200905_wgs_5162_hg38/bwa/result/CF001.chromosome.count"
  
  logger = logging.getLogger('getChromosomeCount')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  print(args)

  count(logger, args.input, args.output)

main()
