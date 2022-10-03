import pysam
import argparse
import sys
import logging
import os
import re
import errno
from asyncore import read

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def write_reads(fout, refmap, reads, counts):
  hostReads = []
  feedReads = []
  hasSecondary = False
  for read in reads:
    if read.is_secondary:
      hasSecondary = True
      continue
    if refmap[read.reference_id]:
      hostReads.append(read)
    else:
      feedReads.append(read)

  hasBoth = False
  savedReads = []
  if len(hostReads) > 0:
    savedReads = hostReads
    if len(feedReads) > 0:
      counts["both"] += 1
    else:
      counts["host"] += 1
  else:
    savedReads = feedReads
    counts["feed"] += 1

  if hasBoth or hasSecondary:
    readCount = {True:0, False:0}
    for read in savedReads:
      readCount[read.is_read1] += 1
    numberOfHit = max(readCount.values())
    #oldNumberOfHit = read.get_tag("NH")
    #print("hasSecondary:%s, hasBoth:%s, old:%d, new:%d" % (str(hasSecondary), str(hasBoth), oldNumberOfHit, numberOfHit))
    for read in savedReads:
      read.set_tag("NH", numberOfHit)

  for read in savedReads:
    fout.write(read)

DEBUG = False
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="filter bam by chromosome priority using pattern.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output BAM file", required=NOT_DEBUG)
parser.add_argument('--host_prefix', action='store', nargs='?', help="Input chromosome prefix for host genome (such as mm10_)", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="/scratch/vickers_lab/projects/20200805_5057_AD_rnaseq_hsammu_combined_byMars.tiger.bam/sort_by_name/result/Blank.sortedByName.bam"
  args.output="/scratch/vickers_lab/projects/20200805_5057_AD_rnaseq_hsammu_combined_byMars.tiger.bam/Blank.name.filtered.bam"
  args.host_prefix="mm10_"
  # args.input="/scratch/vickers_lab/projects/20200805_5057_AD_rnaseq_hsammu_combined_byMars.tiger.bam/chr1.name.bam"
  # args.output="/scratch/vickers_lab/projects/20200805_5057_AD_rnaseq_hsammu_combined_byMars.tiger.bam/chr1.name.filtered.bam"
  # args.host_prefix="mm10_"

logger = logging.getLogger('filterMixBam')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

check_file_exists(args.input)

tmpfile = args.output + ".tmp.bam"

counts={"host":0, "feed":0, "both":0}

refmap = {}
with pysam.Samfile(args.input, "rb") as sam, open(tmpfile, 'w') as output:
  for nf in range(0, sam.nreferences):
    if(sam.get_reference_name(nf) != sam.references[nf]):
      raise Exception("%s != %s" % (sam.get_reference_name(nf), sam.references[nf]))
    refmap[nf] = sam.references[nf].startswith(args.host_prefix)

  header = sam.header
  with pysam.AlignmentFile(tmpfile, "wb", header=header) as fout:
    processed = 0
    accepted = 0
    reads = []
    lastQuery = ""
    for read in sam.fetch(until_eof=True):
      if read.is_unmapped:
        continue

      if read.qname != lastQuery:
        processed += 1
        if processed % 100000 == 0:
          logger.info("processed %d, host %d, feed %d, both %d" % (processed, counts["host"], counts["feed"], counts["both"]))
          
        write_reads(fout, refmap, reads, counts)
        lastQuery = read.qname
        reads = [read]
      else:
        reads.append(read)
    write_reads(fout, refmap, reads, counts)
    logger.info("processed %d, host %d, feed %d, both %d" % (processed, counts["host"], counts["feed"], counts["both"]))

txtfile = os.path.splitext(args.output)[0]+'.txt'
with open(txtfile, "wt") as flog:
  flog.write("Category\tCount\n")
  flog.write("processed\t%d\n" % processed)
  flog.write("host\t%d\n" % counts["host"])
  flog.write("feed\t%d\n" % counts["feed"])
  flog.write("both\t%d\n" % counts["both"])

if os.path.exists(args.output):
  os.unlink(args.output)

os.rename(tmpfile, args.output)
