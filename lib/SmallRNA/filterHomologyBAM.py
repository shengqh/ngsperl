import pysam
import argparse
import sys
import logging
import os
from asyncore import read
from Bio import SeqIO

def outputQuery(outf, saved_read, referencePrefix, homologyPrefix):
  hasRef = any(read.reference_id.startswith(referencePrefix) for read in saved_read)
  if hasRef:
    hasHomo = any(read.reference_id.startswith(homologyPrefix) for read in saved_read)
    if hasRef:
      for read in saved_read:
        if read.reference_id.startswith(homologyPrefix):
          read.reference_id = read.reference_id[length(homologyPrefix):]
          outf.write(read)
  
def filter(outputBAM, inputBAM, referencePrefix, homologyPrefix):
  if inputBAM.endswith(".bam"):
    openmode = "rb"
  else:
    openmode = "r"

  with pysam.AlignmentFile(inputBAM, openmode) as samfile:
    header = sam.header
    with pysam.AlignmentFile(outputBAM, "wb", header=header) as outf:
      processed = 0
      lastQuery =''
      saved_read = []
      for read in samfile.fetch(until_eof=True):
        processed += 1
        if processed % 1000000 == 0:
          logger.info("processed %d" % (processed))
          
        if read.is_unmapped:
          continue;
        
        if lastQuery != read.query_name:
          outputQuery(outf, saved_read, referencePrefix, homologyPrefix)
          saved_read = [read]
          lastQuery = read.query_name
        else:
          saved_read.append(read)
      outputQuery(outf, saved_read, referencePrefix, homologyPrefix)
        
def main():
  parser = argparse.ArgumentParser(description="Filter homology smallRNA mapping BAM file.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = True
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file)', required=NOT_DEBUG)
  parser.add_argument('--referencePrefix', action='store', nargs='?', help='Input reference genome prefix)', required=NOT_DEBUG)
  parser.add_argument('--homologyPrefix', action='store', nargs='?', help='Input homology genome prefix)', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output BAM file", required=NOT_DEBUG)
  
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/cqs/shengq1/temp/RPI11.10000.sam"
    args.referencePrefix="rn5_"
    args.homologyPrefix="mm10_"
    args.output="/scratch/cqs/shengq1/temp/RPI11.rn5.bam"
  
  logger = logging.getLogger('homology')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  filter(args.output, args.input, args.referencePrefix, args.homologyPrefix)
  
if __name__ == "__main__":
    main()
