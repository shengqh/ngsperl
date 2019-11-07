import pysam
import argparse
import sys
import logging
import os
from asyncore import read

parser = argparse.ArgumentParser(description="Parsing BAM file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file', required=True)
parser.add_argument('-g', '--chromosome_size_file', action='store', nargs='?', help='Input chromosome size file', required=True)
parser.add_argument('-f', '--genome_seq_file', action='store', nargs='?', help='Input genome seq file', required=True)
parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output file name", required=True)

args = parser.parse_args()

logger = logging.getLogger('parseBam')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def runCmd(cmd, logger):
  logger.info(cmd)
  os.system(cmd)

bedFile = args.output + ".bed"
dinuFile = args.output + ".dinu.bed"
runCmd("bamToBed -i %s > %s" % (args.input, bedFile), logger)
runCmd("bedtools flank -i %s -g %s -l 2 -r 0 -s > %s" % (bedFile, args.chromosome_size_file, dinuFile), logger)
runCmd("bedtools getfasta -fi %s -bed %s -s -fo %s" % (args.genome_seq_file, dinuFile, args.output), logger)
