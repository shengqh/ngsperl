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

parser = argparse.ArgumentParser(description="Preprocess rose2",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input peak bed file', required=NotDEBUG)
parser.add_argument('-r', '--reads', action='store', nargs='?', help='Input read bam files', required=NotDEBUG)
parser.add_argument('-c', '--control', action='store', nargs='?', help='Input control bam files', required=NotDEBUG)
parser.add_argument('-g', '--genome', action='store', nargs='?', help='Input genome (HG38/HG19/MM10)', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output prefix", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/jbrown_lab/shengq2/projects/20201110_chipmentation_SampleID_tfgene/0_uncompress_peak/result/GYZ241_NAFLD_H3K27ac_rep1.peak.bed"
  args.reads = "/scratch/data/vumc_h3K27ac_pilot/chipmentation/DRX0FMN_liverPilot_200701/4615-LD-1/align/rep1/4615-LD-1_S57_L005_R1_001.nodup.bam,/scratch/data/vumc_h3K27ac_pilot/chipmentation/DRX0FMN_liverPilot_200701/4615-LD-2/align/rep1/4615-LD-2_S58_L005_R1_001.nodup.bam"
  args.control = "/scratch/data/vumc_h3K27ac_pilot/chipmentation/DRX0FMN_liverPilot_200701/4615-LD-1/align/ctl1/4615-LD-3_S59_L005_R1_001.nodup.bam,/scratch/data/vumc_h3K27ac_pilot/chipmentation/DRX0FMN_liverPilot_200701/4615-LD-2/align/ctl1/4615-LD-3_S59_L005_R1_001.nodup.bam"
  args.genome = "HG38"
  args.output = "GYZ241_NAFLD_H3K27ac_rep1/"

logger = logging.getLogger('rose2')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

reads_files = args.reads.split(',')
cmd = "ROSE2 -g %s -i %s -r %s" % (args.genome, args.input, reads_files[0])
if len(reads_files) > 1:
  cmd = cmd + " -b %s" % ",".join(reads_files[1:])

control_files = args.control.split(',')
cmd = cmd + " -c %s" % control_files[0]
cmd = cmd + " -o %s" % args.output

runCmd(cmd, logger)

logger.info("done.")