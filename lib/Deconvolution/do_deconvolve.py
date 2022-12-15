import argparse
import sys
import logging
import os
import os.path
import logging
import errno
import shutil
from collections import OrderedDict
import deconvolve as deconv

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

parser = argparse.ArgumentParser(description="deconvolve",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input normalized cpm csv file', required=NotDEBUG)
parser.add_argument('-b', '--basic', action='store', nargs='?', help='Input basic matrix file (tsp_v1_basisMatrix.txt)', required=NotDEBUG)
parser.add_argument('-n', '--name', action='store', nargs='?', help='Input sample name in cpm file', required=NotDEBUG)

args = parser.parse_args()

logger = logging.getLogger('deconvolve')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

deconv.main(1,  args.input, args.basic, [args.name] , 'nuSVR', args.name)