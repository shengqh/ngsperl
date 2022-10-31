import argparse
import logging
import os
import os.path
import logging
import errno
import re;
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

parser = argparse.ArgumentParser(description="Get count table",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count file list', required=NotDEBUG)
parser.add_argument('-p', '--pattern', action='store', nargs='?', help='Input feature pattern', required=NotDEBUG)
parser.add_argument('--noheader', action='store_true', help='If the count file without header')
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NotDEBUG)

args = parser.parse_args()

logger = logging.getLogger('count_table')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

count_files = readFileMap(args.input)

count_table = OrderedDict()

samples = list(count_files.keys())

for sample in count_files.keys():
  count_file = count_files[sample]
  logger.info(f"reading {count_file}")

  with open(count_file, "rt") as fin:
    if not args.noheader:
      fin.readline()
    for line in fin:
      name, count = line.rstrip().split('\t')
      name = name.replace('"', '')
      if re.match(args.pattern, name):
        if not name in count_table:
          count_table[name] = OrderedDict()
        count_table[name][sample] = count

with open(args.output, "wt") as fout:
  sample_str = "\t".join(samples)
  fout.write(f"Feature\t{sample_str}\n")
  for feature in count_table.keys():
    fout.write(feature)
    sample_count = count_table[feature]
    for sample in samples:
      if sample in sample_count:
        fout.write(f"\t{sample_count[sample]}")
      else:
        fout.write(f"\t0")
    fout.write("\n")

logger.info("done.")