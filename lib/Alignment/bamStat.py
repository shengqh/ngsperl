import pysam
import argparse
import sys
import logging
import os
import re
from asyncore import read

from enum import Enum

class ReadType(Enum):
  UNMAPPED = 1
  UNIQUELY_MAPPED = 2
  MULTIPLE_MAPPED = 3

def get_read_type(reads):
  read1_count = 0
  read2_count = 0

  for read in reads:
    if read.is_unmapped:
      return (ReadType.UNMAPPED)

    if read.is_read1:
      read1_count += 1
    else:
      read2_count += 1
  
  if (read1_count < 2) and (read2_count < 2):
    return(ReadType.UNIQUELY_MAPPED)
  else:
    return(ReadType.MULTIPLE_MAPPED)

DEBUG = False
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="Get bam summary information",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file, sorted by query name', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output summary file", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="/scratch/vickers_lab/projects/20200805_5057_AD_rnaseq_hsammu_combined_byMars.tiger.bam/sort_by_name/result/Blank.sortedByName.bam"
  args.output=args.input + ".stat"

logger = logging.getLogger('bamSummary')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

with pysam.Samfile(args.input, "rb") as sam:
  type_dic = {ReadType.UNMAPPED:0, ReadType.UNIQUELY_MAPPED:0, ReadType.MULTIPLE_MAPPED:0}

  read_types = {}

  processed = 0
  last_query = ""
  reads = []
  for read in sam.fetch(until_eof=True):
    processed += 1
    if processed % 1000000 == 0:
      logger.info("reads %d = fragments %d: unmapped %d, uniquely mapped %d, multiple_mapped %d" % 
        (processed, sum(type_dic.values()), type_dic[ReadType.UNMAPPED], type_dic[ReadType.UNIQUELY_MAPPED], type_dic[ReadType.MULTIPLE_MAPPED] ))

    if read.qname == last_query:
      reads.append(read)
      continue

    if len(reads) > 0:
      type_dic[get_read_type(reads)] += 1

    reads = [read]
    last_query = read.qname

#last read group
type_dic[get_read_type(reads)] += 1

with open(args.output, "wt") as fout:
  fout.write("Category\tCount\n")
  fout.write("TotalEntries\t%d\nTotalFragments\t%d\nUnmappedFragments\t%d\nUniquelyMappedFragments\t%d\nMultipleMappedFragments\t%d\n" % 
        (processed, sum(type_dic.values()), type_dic[ReadType.UNMAPPED], type_dic[ReadType.UNIQUELY_MAPPED], type_dic[ReadType.MULTIPLE_MAPPED] ))
  
logger.info("done")