#!/usr/bin/env python3

import os
import os.path
import argparse
import sys
import logging
import re
 
from glob import glob
from clean_utils import *

def clean_mutect(logger, project_dir, remove_patterns, remove=True):
  log_folder =  os.path.join(project_dir, "bwa_refine_muTect_01_call", "log")
  samples = find_error_samples_by_mutect(logger, log_folder)
  #print(samples)
  remove_files(logger, project_dir, samples, remove_patterns, remove=remove)
  return(samples)

def clean_mutect2(logger, project_dir, remove_patterns, remove=True):
  log_folder =  os.path.join(project_dir, "bwa_refine_muTect2indel", "log")
  samples = find_error_samples_by_mutect2(logger, log_folder)
  #print(samples)
  remove_files(logger, project_dir, samples, remove_patterns, remove=remove)
  return(samples)

def main():
  parser = argparse.ArgumentParser(description="Clean RNAseq project",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = True
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input project root folder', required=NOT_DEBUG)
  parser.add_argument('-d', '--delete', action='store_true', help='Delete files')
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    #args.input="/scratch/stein_lab/20210429_qiong_rnaseq_6130_mir1246_transfection_hg38/"
    args.input="/scratch/cqs/PCA_scRNAseq/Exoseq/20210430_6109_YX"

  args.input= os.path.abspath(args.input)
  print(args)
  
  logger = logging.getLogger('clean_mutect2')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  remove_mutect_files = [
    ["bwa_refine_muTect_01_call/result/*", ".somatic.*"],
    ["bwa_refine_muTect_01_call/log", "_mt.log"],
  ]

  remove_mutect2_files = [
    ["bwa_refine_muTect2indel/result", ".somatic.*"],
    ["bwa_refine_muTect2indel/log", "_mt2.log"],
  ]

  clean_mutect(logger, project_dir=args.input, remove_patterns=remove_mutect_files, remove=args.delete)
  clean_mutect2(logger, project_dir=args.input, remove_patterns=remove_mutect2_files, remove=args.delete)

if __name__ == "__main__":
    main()
