#!/usr/bin/env python3

import os
import os.path
import argparse
import sys
import logging
import re
 
from glob import glob
from clean_utils import *

def check_star(logger, project_dir, remove_patterns, delete):
  star_folder =  os.path.join(project_dir, "star_featurecount.v2")
  samples = find_error_samples_by_star(logger, star_folder)
  if delete:
    remove_files(logger, project_dir, samples, remove_patterns)
  else:
    print(samples)
  
  return(samples)

def main():
  parser = argparse.ArgumentParser(description="Clean RNAseq project",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input project root folder', required=NOT_DEBUG)
  parser.add_argument('-d', '--delete', action='store_true', help='Delete files')
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    #args.input="/scratch/stein_lab/20210429_qiong_rnaseq_6130_mir1246_transfection_hg38/"
    args.input="/scratch/stein_lab/20210429_qiong_rnaseq_6130_gut_hg38"

  args.input= os.path.abspath(args.input)
  print(args)
  
  logger = logging.getLogger('check_star')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  remove_star_files = [
    ["star_featurecount.v2/result", ".count*"],
    ["star_featurecount.v2/result", ".bamstat"],
    ["star_featurecount.v2/result", ".00??.bam"],
    ["star_featurecount.v2/result", "__*"],
    ["star_featurecount.v2/result", "_Aligned*"],
    ["star_featurecount.v2/result", "_Log*"],
    ["star_featurecount.v2/result", "_SJ.out.tab"],
    ["star_featurecount.v2/result", ".splicing.bed"],
    ["star_featurecount.v2/log", "_sf.log"],
  ]

  #clean_by_posttrim_fastqc(logger, project_dir=args.input, remove_patterns=remove_posttrim_fastqc_files)
  #clean_by_posttrim_fastqc_and_fastq_len(logger, project_dir=args.input, remove_patterns=remove_cutadapt_files+remove_posttrim_fastqc_files+remove_star_files+remove_posttrim_fastq_len_files)

  #clean_by_cutadapt(logger, project_dir=args.input, remove_patterns=remove_cutadapt_files)
  #clean_by_posttrim_fastqc(logger, project_dir=args.input, remove_patterns=remove_cutadapt_files+remove_posttrim_fastqc_files+remove_star_files+remove_posttrim_fastq_len_files)
  #clean_fastq_len(logger, project_dir=args.input, remove_patterns=remove_cutadapt_files+remove_posttrim_fastqc_files+remove_star_files+remove_posttrim_fastq_len_files)
  check_star(logger, project_dir=args.input, remove_patterns=remove_star_files, delete=args.delete)

if __name__ == "__main__":
    main()
