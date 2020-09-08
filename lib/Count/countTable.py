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
      check_file_exists(filepath.strip())
  return(result)

DEBUG=True
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Get count table",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count file list', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/jbrown_lab/shengq2/projects/20200905_wgs_5162_hg38/bwa_summary/result/wgs_5162__fileList1.list"
  args.output = "/scratch/jbrown_lab/shengq2/projects/20200905_wgs_5162_hg38/bwa_summary/result/wgs_5162.txt"

logger = logging.getLogger('countTable')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

count_files = readFileMap(args.input)
counts = OrderedDict()

for name in count_files.keys():
  count_file = count_files[name]
  with open(count_file, "rt") as fin:
    count = 0
    for line in fin:
      if line.startswith('#'):
        continue
      parts = line.rstrip().split('\t')
      if parts[1].isnumeric():
        counts.setdefault(parts[0], {})[name] = parts[1]

with open(args.output, "wt") as fout:
  fout.write("Feature\t%s\n" % ("\t".join(count_files.keys())))
  for feature in counts.keys():
    fout.write("%s\t%s\n" % (feature, "\t".join([counts[feature][name] for name in count_files.keys()])))

logger.info("done.")