import pysam
import argparse
import sys
import logging
import os
from asyncore import read
from Bio import SeqIO

def outputQuery(outf, saved_read, referencePrefix, chrMap):
  for read in saved_read:
    if chrMap[read.reference_id].startswith(referencePrefix):
      read.reference_id = chrMap[chrMap[read.reference_id]]
      outf.write(read)
  
def filter(outputBAM, inputBAM, referencePrefix, logger):
  if inputBAM.endswith(".bam"):
    openmode = "rb"
  else:
    openmode = "r"

  with pysam.AlignmentFile(inputBAM, openmode) as samfile:
    header = samfile.header
    chrMap = {}
    sqs = header["SQ"]
    for idx in range(len(sqs)-1, -1, -1):
      sq = sqs[idx]
      chr = sq['SN']
      chrMap[idx] = chr
      if not chr.startswith(referencePrefix):
        del sqs[idx]
    
    for idx in range(0, len(sqs)):
      sq = sqs[idx]
      chr = sq['SN']
      chrMap[chr] = idx
      sq['SN'] = chr[len(referencePrefix):]
      
    header["SQ"] = sqs
      
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
          outputQuery(outf, saved_read, referencePrefix, chrMap)
          saved_read = [read]
          lastQuery = read.query_name
        else:
          saved_read.append(read)
      outputQuery(outf, saved_read, referencePrefix, chrMap)
        
def main():
  parser = argparse.ArgumentParser(description="Filter reference smallRNA mapping BAM file.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file', required=NOT_DEBUG)
  parser.add_argument('--referencePrefix', action='store', nargs='?', help='Input reference genome prefix', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output BAM file", required=NOT_DEBUG)
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/cqs/shengq1/temp/RPI11.10000.sam"
    args.referencePrefix="rn5_"
    args.output="/scratch/cqs/shengq1/temp/RPI11.10000.rn5_.bam"
  
  logger = logging.getLogger('reference')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  filter(args.output, args.input, args.referencePrefix, logger)
  
if __name__ == "__main__":
    main()
