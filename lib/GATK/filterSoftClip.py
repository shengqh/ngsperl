import pysam
import argparse
import sys
import logging
import os
from asyncore import read

parser = argparse.ArgumentParser(description="Remove low quality or soft-clip reads from bam file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file', required=True)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output bam file", required=True)
parser.add_argument('--min-mapq', action='store', nargs='?', type=int, default=10, help="Minimum mapping quality of read")
#parser.add_argument('-m', '--minimum_softclip_bases', action='store', nargs='?', default=3, help="Discard reads with minimum softclip bases")

args = parser.parse_args()

logger = logging.getLogger('filterSoftClip')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

if args.input.endswith(".bam"):
  openmode = "rb"
else:
  openmode = "r"

def filterReadQuality(read, min_mapq):
  return(read.is_unmapped or read.mapq < min_mapq or read.is_secondary or read.is_qcfail or read.is_duplicate or read.is_supplementary)

def filterSoftClip(read):
  return("S" in read.cigarstring)

def filterSoftClipNumber(read):
  #https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment
  ctf = read.cigartuples[0]
  if ctf[0] == 4: #soft clip
    if ctf[1] >= args.minimum_softclip_bases:
      return(True)
  
  if len(read.cigartuples) > 1:
    ctf = read.cigartuples[len(read.cigartuples)-1]
    if ctf[0] == 4: #soft clip
      if ctf[1] >= args.minimum_softclip_bases:
        return(True)
  
  return(False)

logger.info("Finding reads with soft clip ...")
filterQuerySet = set()
with pysam.AlignmentFile(args.input, openmode) as samfile:
  processed = 0
  for read in samfile.fetch(until_eof=True):
    processed += 1
    if processed % 1000000 == 0:
      logger.info("processed %d entries, found %d queries with soft clip ..." % (processed, len(filterQuerySet)))

    if read.is_unmapped:
      continue;
      
    if filterSoftClip(read):
      filterQuerySet.add(read.qname)
      continue

logger.info("Found %d queries with soft clip ..." % len(filterQuerySet))
        
tmpfile = args.output + ".tmp.bam"
logger.info("Filtering reads with low mapping quality or soft clip ...")
with pysam.AlignmentFile(args.input, openmode) as samfile:
  with pysam.AlignmentFile(tmpfile, "wb", header=samfile.header) as outf:
    processed = 0
    written = 0
    for read in samfile.fetch(until_eof=True):
      processed += 1
      if processed % 1000000 == 0:
        logger.info("processed %d entries, written %d entries ..." % (processed, written))

      if filterReadQuality(read, args.min_mapq):
        continue;

      if read.qname in filterQuerySet:
        continue

      written += 1
      outf.write(read)

    logger.info("totally processed %d, written %d" % (processed, written))

if os.path.isfile(args.output):
  os.remove(args.output)
os.rename(tmpfile, args.output)

logger.info("Done")

