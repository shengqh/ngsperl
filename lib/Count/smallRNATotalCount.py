import argparse
import sys
import logging
import os
import os.path
import logging
import errno
import shutil
import pysam

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Get total count",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bam file list', required=NotDEBUG)
parser.add_argument('-c', '--count', action='store', nargs='?', help='Input count file list', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/vickers_lab/projects/20200924_tellseq/pseudohap_bowtie1_hg38/result/HDL_96.bam"
  args.count = "/scratch/vickers_lab/projects/20200902_5188_ES_smRNA_human_v5_byMars/preprocessing/identical/result/HDL_96_clipped_identical.fastq.dupcount"
  args.output = "HDL_96.pseudohap.hg38.count"

logger = logging.getLogger('totalCount')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

check_file_exists(args.input)
check_file_exists(args.count)

logger.info(f"reading {args.count} ...")
countMap = {}
with open(args.count, "rt") as fin:
  header = fin.readline()
  for line in fin:
    parts = line.split('\t')
    countMap[parts[0]] = int(parts[1])

qnames = set()
logger.info(f"processing {args.input} ...")
totalCount = 0
totalPMCount = 0
with pysam.AlignmentFile(args.input, "rb") as sf:
  for s in sf.fetch():
    if s.is_unmapped:
      continue
    if s.query_name in qnames:
      continue
    nm = s.get_tag("NM")
    curCount = countMap[s.query_name]
    totalCount += curCount

    if nm == 0:
      totalPMCount += curCount
      
    qnames.add(s.query_name)

with open(args.output, "wt") as fout:
  fout.write(f"Total\t{totalCount}\n")
  fout.write(f"TotalPM\t{totalPMCount}\n")

logger.info("done.")