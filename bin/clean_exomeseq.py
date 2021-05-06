#!/usr/bin/env python3

import os
import os.path
import argparse
import sys
import logging
import re
 
from glob import glob
from clean_utils import *

def clean_bwa(logger, project_dir, remove_patterns, remove=True):
  bwa_folder =  os.path.join(project_dir, "bwa")
  samples = find_error_samples_by_bwa(logger, bwa_folder)
  #print(samples)
  remove_files(logger, project_dir, samples, remove_patterns, remove=remove)
  return(samples)

def main():
  parser = argparse.ArgumentParser(description="Clean RNAseq project",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input project root folder', required=NOT_DEBUG)
  parser.add_argument('-l', '--list', action='store_true', help='List all files will be deleted')
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    #args.input="/scratch/stein_lab/20210429_qiong_rnaseq_6130_mir1246_transfection_hg38/"
    args.input="/scratch/stein_lab/20210429_qiong_rnaseq_6130_gut_hg38"

  args.input= os.path.abspath(args.input)
  print(args)
  
  logger = logging.getLogger('clean_rnaseq')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  remove_bwa_files = [
    ["bwa/result", ".bamstat"],
    ["bwa/result", ".bwa.version"],
    ["bwa/result", ".sortedByCoord.bam*"],
    ["bwa/log", "_bwa.log"],
  ]

  remove_bwa_refine_files = [
    ["bwa_refine/result", ".rmdup.*"],
    ["bwa_refine/log", "_rf.log"],
  ]

  clean_bwa(logger, project_dir=args.input, remove_patterns=remove_bwa_files+remove_bwa_refine_files, remove=not args.list)

if __name__ == "__main__":
    main()
