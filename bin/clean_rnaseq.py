#!/usr/bin/env python3

import os
import os.path
import argparse
import sys
import logging
import re
 
from glob import glob
from clean_utils import *

def clean_by_cutadapt(logger, project_dir, remove_patterns):
  cutadapt_result_folder = os.path.join(project_dir, "cutadapt", "result")
  samples = do_clean_by_cutadapt(logger, project_dir, cutadapt_result_folder, remove_patterns)
  return(samples)

def clean_by_posttrim_fastqc(logger, project_dir, remove_patterns):
  post_trim_fastqc_log_folder = os.path.join(project_dir, "fastqc_post_trim", "log")
  samples = find_error_samples_by_fastqc(logger, post_trim_fastqc_log_folder)
  remove_files(logger, project_dir, samples, remove_patterns)
  return(samples)

def clean_fastq_len(logger, project_dir, remove_patterns):
  fastq_len_result_folder = os.path.join(project_dir, "fastq_len", "result")
  samples = find_error_samples_by_fastq_len(logger, fastq_len_result_folder)
  remove_files(logger, project_dir, samples, remove_patterns)
  return(samples)

def clean_star(logger, project_dir, remove_patterns):
  star_folder =  os.path.join(project_dir, "star_featurecount")
  samples = find_error_samples_by_star(logger, star_folder)
  #print(samples)
  remove_files(logger, project_dir, samples, remove_patterns)
  return(samples)

def main():
  parser = argparse.ArgumentParser(description="Clean RNAseq project",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input project root folder', required=NOT_DEBUG)
  
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
  
  remove_cutadapt_files = [
    ["cutadapt/result", "_clipped*"],
    ["cutadapt/result", ".version"],
    ["cutadapt/log", "_cut.log"],
  ]

  remove_posttrim_fastqc_files = [
    ["fastqc_post_trim/result", "/*"],
    ["fastqc_post_trim/log", "_fq.log"],
  ]

  remove_posttrim_fastq_len_files = [
    ["fastq_len/result", ".len"],
    ["fastq_len/result", ".len.error"],
    ["fastq_len/log", "_flen.log"],
  ]

  remove_star_files = [
    ["star_featurecount/result", ".count*"],
    ["star_featurecount/result", ".bamstat"],
    ["star_featurecount/result", "__*"],
    ["star_featurecount/result", "_Aligned*"],
    ["star_featurecount/result", "_Log*"],
    ["star_featurecount/result", "_SJ.out.tab"],
    ["star_featurecount/result", ".splicing.bed"],
    ["star_featurecount/log", "_sf.log"],
  ]

  clean_by_cutadapt(logger, project_dir=args.input, remove_patterns=remove_cutadapt_files)
  clean_by_posttrim_fastqc(logger, project_dir=args.input, remove_patterns=remove_cutadapt_files+remove_posttrim_fastqc_files+remove_star_files+remove_posttrim_fastq_len_files)
  clean_fastq_len(logger, project_dir=args.input, remove_patterns=remove_cutadapt_files+remove_posttrim_fastqc_files+remove_star_files+remove_posttrim_fastq_len_files)
  clean_star(logger, project_dir=args.input, remove_patterns=remove_star_files)

if __name__ == "__main__":
    main()
