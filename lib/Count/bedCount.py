import argparse
import sys
import logging
import os
import os.path
import logging
import errno
import shutil
from collections import OrderedDict

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def readFileMap(fileName):
  check_file_exists(fileName)

  result = OrderedDict()
  with open(fileName) as fh:
    for line in fh:
      filepath, name = line.strip().split('\t', 1)
      result[name] = filepath.strip()
  return(result)

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Get count in bed files",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input peak bed file list', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/jbrown_lab/shengq2/projects/20200720_chipseq_GSE140641_mouse/croo/result/H3K27ac_activated/peak/overlap_reproducibility/overlap.optimal_peak.narrowPeak.gz.bed"
  args.output = "H3K27ac_activated.txt"

logger = logging.getLogger('bedCount')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

bed_files = readFileMap(args.input)
bed_counts = OrderedDict()

for bed_name in bed_files.keys():
  bed_file = bed_files[bed_name]
  with open(bed_file, "rt") as fin:
    count = 0
    for line in fin:
      if not line.startswith('#'):
        count = count + 1
    bed_counts[bed_name] = count

with open(args.output, "wt") as fout:
  fout.write("Sample\tPeakCount\n")
  for bed_name in bed_counts.keys():
    fout.write("%s\t%d\n" % (bed_name, bed_counts[bed_name]))

logger.info("done.")