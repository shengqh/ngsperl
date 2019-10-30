import pysam
import argparse
import sys
import logging
import os
from asyncore import read

parser = argparse.ArgumentParser(description="Build soft clip position distribution in BAM file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

DEBUG=False
NOT_DEBUG = not DEBUG

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file', required=NOT_DEBUG)
parser.add_argument('--min-mapq', action='store', nargs='?', type=int, default=10, help="Minimum mapping quality of read")
parser.add_argument('--binsize', action='store', nargs='?', type=int, default=1000, help="Bin size of position")
parser.add_argument('--min-depth', action='store', nargs='?', type=int, default=100, help="Minimum depth for output")
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output soft clip distribution file name", required=NOT_DEBUG)

if NOT_DEBUG and len(sys.argv)==1:
  parser.print_help()
  sys.exit(1)
  
args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/cqs/shengq2/jennifer/20190906_lindsay_exomeseq_3772_hg38/softclip/P_175_06.indel.recal.TP53.bam"
  args.output = "/scratch/cqs/shengq2/jennifer/20190906_lindsay_exomeseq_3772_hg38/softclip/P_175_06.softclip.position.tsv"

logger = logging.getLogger('bamSoftClipPosition')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def filterReadQuality(read, min_mapq):
  return(read.is_unmapped or read.mapping_quality < min_mapq or read.is_secondary or read.is_qcfail or read.is_duplicate or read.is_supplementary)

def hasSoftClip(read):
  return("S" in read.cigarstring)

chrPositionMap = {}
processed = 0
logger.info("reading %s" % args.input)
with pysam.Samfile(args.input, "rb") as samfile:
  for read in samfile.fetch(until_eof=True):
    processed += 1
    if processed % 1000000 == 0:
      logger.info("processed %d" % processed)
      #break
      
    if filterReadQuality(read, args.min_mapq):
      continue
    
    if len(read.reference_name) > 5:
      continue

    if not read.reference_name in chrPositionMap:
      chrPositionMap[read.reference_name] = {}
    positionMap = chrPositionMap[read.reference_name]

    position = int(read.reference_start / args.binsize)
    if not position in positionMap:
      positionMap[position] = [0, 0]
    posvalues = positionMap[position]

    if hasSoftClip(read):
      posvalues[0] = posvalues[0] + 1
    else:
      posvalues[1] = posvalues[1] + 1

with open(args.output, "wt") as sw:  
  sw.write("Chr\tStartPosition\tSoftClipRead\tOtherRead\tSoftClipPerc\n")
  for chr in chrPositionMap.keys():
    positionMap = chrPositionMap[chr]
    positions = sorted(positionMap.keys())
    for pos in positions:
      posvalues = positionMap[pos]
      sread = posvalues[0]
      oread = posvalues[1]
      allread = sread + oread
      if allread >= args.min_depth:
        sw.write("%s\t%d\t%d\t%d\t%.2f\n" % (chr, pos * args.binsize, sread, oread, sread * 1.0 / allread) )

logger.info("done.")