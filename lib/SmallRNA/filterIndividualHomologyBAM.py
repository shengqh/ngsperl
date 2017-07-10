import pysam
import argparse
import sys
import logging
import os
from asyncore import read
from Bio import SeqIO

def getNumberOfMismatch(read):
  result = 0;
  if read.has_tag("NM"):
    result = read.get_tag("NM")
  elif read.has_tag("nM"):   
    result = read.get_tag("nM")
  return(result)
  
def filter(outputBAM, inputReferenceBAM, inputHomologyBAM, logger):
  if inputReferenceBAM.endswith(".bam"):
    openmode = "rb"
  else:
    openmode = "r"

  #get query and number of mismatch map in reference BAM
  refQuery = {}      
  with pysam.AlignmentFile(inputReferenceBAM, openmode) as samfile:
    processed = 0
    for read in samfile.fetch(until_eof=True):
      processed += 1
      if processed % 1000000 == 0:
        logger.info("read reference entries %d" % (processed))
          
      if read.is_unmapped:
        continue;
      
      refQuery[read.query_name] = getNumberOfMismatch(read);

  #filter homology BAM
  if inputHomologyBAM.endswith(".bam"):
    openmode = "rb"
  else:
    openmode = "r"

  with pysam.AlignmentFile(inputHomologyBAM, openmode) as samfile:
    header = samfile.header
    with pysam.AlignmentFile(outputBAM, "wb", header=header) as outf:
      processed = 0
      for read in samfile.fetch(until_eof=True):
        processed += 1
        if processed % 1000000 == 0:
          logger.info("processed %d" % (processed))
          
        if read.is_unmapped:
          continue;
        
        if not read.query_name in refQuery:
          continue;
          
        refNumberOfMismatch = refQuery[read.query_name];
        homoNumberOfMismatch = getNumberOfMismatch(read);
        
        if(refNumberOfMismatch > homoNumberOfMismatch):
          continue;
          
        outf.write(read);
        
def main():
  parser = argparse.ArgumentParser(description="Filter homology smallRNA mapping BAM file.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('--referenceBAM', action='store', nargs='?', help='Input reference BAM', required=NOT_DEBUG)
  parser.add_argument('--homologyBAM', action='store', nargs='?', help='Input homology BAM', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output filtered homology BAM", required=NOT_DEBUG)
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    args.referenceBAM="/scratch/cqs/shengq1/vickers/20170628_smallRNA_2795_2516_rat_homology/star_rn5/result/Rat_HDL_Lean1_12_Aligned.out.bam"
    args.homologyBAM="/scratch/cqs/shengq1/vickers/20170628_smallRNA_2795_2516_rat_homology/star_mm10/result/Rat_HDL_Lean1_12_Aligned.out.bam"
    args.output="/scratch/cqs/shengq1/vickers/20170628_smallRNA_2795_2516_rat_homology/star_mm10_filtered/result/Rat_HDL_Lean1_12.filtered.bam"
  
  logger = logging.getLogger('homology')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  filter(args.output, args.referenceBAM, args.homologyBAM, logger)
  
if __name__ == "__main__":
    main()
