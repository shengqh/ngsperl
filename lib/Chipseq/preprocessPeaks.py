import argparse
import sys
import logging
import os
import pybedtools 
import requests
import shutil

def runCmd(cmd, logger):
  logger.info(cmd)
  os.system(cmd)

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Preprocess peaks",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input gzipped peak bed file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output prefix", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/data/vumc_h3K27ac_pilot/chipmentation/DRX0FMN_liverPilot_200701/4615-LD-1/peak/overlap_reproducibility/overlap.conservative_peak.narrowPeak.gz,/scratch/data/vumc_h3K27ac_pilot/chipmentation/DRX0FMN_liverPilot_200701/4615-LD-2/peak/overlap_reproducibility/overlap.conservative_peak.narrowPeak.gz"
  args.output = "GYZ241_NAFLD_H3K27ac_rep1.peaks.bed"

logger = logging.getLogger('preprocessPeaks')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

peaks_files = args.input.split(',')
if len(peaks_files) > 1:
  otmpfile = args.output + ".tmp"
  unzip_files = []
  for i in range(len(peaks_files)):
    pf = peaks_files[i]
    tmpfile = "%s_%d.bed" % (args.output, i)
    cmd = "gunzip -c %s > %s" % (pf, tmpfile)
    runCmd(cmd, logger)
    unzip_files.append(tmpfile)
  cmd = "bedops --merge %s > %s" % (" ".join(unzip_files), otmpfile)
  runCmd(cmd, logger)
  cmd = "rm %s" % " ".join(unzip_files)
  runCmd(cmd, logger)

  with open(otmpfile, "rt") as fin, open(args.output, "wt") as fout:
    index = 0
    for line in fin:
      index += 1
      fout.write("%s\tPeak_%d\n" % (line.rstrip(), index))

  os.remove(otmpfile)
else:
  cmd = "gunzip -c %s > %s" % (args.input, args.output)
  runCmd(cmd, logger)

logger.info("done.")