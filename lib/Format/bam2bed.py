import pysam
import argparse
import sys
import logging
import os
from asyncore import read

#this bam2bed is specific designed for DNASeq without structure variation
#filter paired-end data by criteria used in MACS2, https://groups.google.com/forum/#!topic/macs-announcement/Vh916L3tOOE
parser = argparse.ArgumentParser(description="Convert SAM/BAM to bed. \nFor paired end data, only the fragment length between min-fragment-length and max-fragment-length will be exported. The reference start will be the minimum start of read1 and read2 and the reference end will be reference start plus fragment length. The strand will be set as same as read1.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file (use "-" as stdin)', required=True)
parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output bed file name")
parser.add_argument('--min-mapq', action='store', nargs='?', type=int, default=10, help="Minimum mapping quality of read")
parser.add_argument('--min-fragment-length', action='store', nargs='?', type=int, default=30, help="Minimum fragment length of paired end reads")
parser.add_argument('--max-fragment-length', action='store', nargs='?', type=int, default=1000, help="Maximum fragment length of paired end reads")
parser.add_argument('--shift_forward', action='store', nargs='?', type=int, default=0, help="Shift bases for forward strand read")
parser.add_argument('--shift_reverse', action='store', nargs='?', type=int, default=0, help="Shift bases for reverse strand read")

args = parser.parse_args()

logger = logging.getLogger('bam2bed')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

if args.input.endswith(".bam"):
  openmode = "rb"
else:
  openmode = "r"

if args.output == "-":
  output = sys.stdout
else:
  tmpfile = args.output + ".tmp"
  output = open(tmpfile, 'w')


samfile = pysam.Samfile(args.input, openmode)
try:
  processed = 0
  accepted = 0
  saved_read = {}
  for read in samfile.fetch(until_eof=True):
    processed += 1
    if processed % 1000000 == 0:
      logger.info("processed %d, accepted %d, cache %d" % (processed, accepted, len(saved_read)))
      
    if read.is_unmapped or read.mapq < args.min_mapq or read.is_secondary or read.is_qcfail or read.is_duplicate or read.is_supplementary:
      continue;
    
    if read.is_paired:
      if not read.is_proper_pair or read.mate_is_unmapped:
        continue;
      
      if read.is_reverse == read.mate_is_reverse:
        continue;
      
      if read.reference_name != read.next_reference_name:
        continue;
      
      if not saved_read.has_key(read.qname):
        saved_read[read.qname] = read
        continue
      
      paired_read = saved_read[read.qname]
      
      assert read.is_read1 != paired_read.is_read1, "Something went wrong (two records indicates same read1 or read2: %s), ignored" % read.qname
      
      try:
        reference_start = min(read.reference_start, paired_read.reference_start)
        reference_end = max(read.reference_end, paired_read.reference_end)
        
        fragment_length = reference_end - reference_start
        if fragment_length < args.min_fragment_length or fragment_length > args.max_fragment_length:
          continue
        
        mapq = min(read.mapq, paired_read.mapq)
        is_reverse = (read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse)
      finally:
        del saved_read[read.qname]
    else:
      is_reverse = read.is_reverse
      reference_start = read.reference_start
      reference_end = read.reference_end
      mapq = read.mapq
  
    if is_reverse:    
      strand = "-"
      shift = args.shift_reverse
    else:
      strand = "+"
      shift = args.shift_forward
    output.write("%s\t%d\t%d\t%s\t%d\t%s\n" % (read.reference_name, reference_start + shift, reference_end + shift, read.qname, mapq, strand))
    output.flush()
    accepted += 1
  
  logger.info("totally processed %d, accepted %d, ignore_unpaired %d" % (processed, accepted, len(saved_read)))
  samfile.close()
  output.close()
  if args.output != "-":
    if os.path.isfile(args.output):
      os.delete(args.output)
    os.rename(tmpfile, args.output)
except:
  samfile.close()
  output.close()
  raise
        
  
