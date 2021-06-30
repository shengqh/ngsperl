#!/usr/bin/env python3

import os
import os.path
import argparse
import sys
import logging
import re
 
from glob import glob
from clean_utils import *


def clean_by_cutadapt(logger, project_dir, remove_patterns, remove=True):
  cutadapt_result_folder = os.path.join(project_dir, "cutadapt", "result")
  samples = do_clean_by_cutadapt(logger, project_dir, cutadapt_result_folder, remove_patterns, remove=remove)
  return(samples)

def clean_by_posttrim_fastqc(logger, project_dir, remove_patterns, remove=True):
  post_trim_fastqc_log_folder = os.path.join(project_dir, "fastqc_post_trim", "log")
  samples = find_error_samples_by_fastqc(logger, post_trim_fastqc_log_folder)
  remove_files(logger, project_dir, samples, remove_patterns, remove=remove)
  return(samples)

def clean_bwa(logger, project_dir, remove_patterns, remove=True):
  bwa_folder =  os.path.join(project_dir, "bwa")
  samples = find_error_samples_by_bwa(logger, bwa_folder)
  #print(samples)
  remove_files(logger, project_dir, samples, remove_patterns, remove=remove)
  return(samples)

def clean_bwa_refine(logger, project_dir, remove_patterns, remove=True):
  folder =  os.path.join(project_dir, "bwa_refine")
  samples = find_error_samples_by_bwa_refine(logger, folder)
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

  remove_bwa_files = [
    ["bwa/result", ".bamstat"],
    ["bwa/result", ".bwa.version"],
    ["bwa/result", ".sortedByCoord.bam*"],
    ["bwa/result", ".unsorted.bam*"],
    ["bwa/result", ".bwa.stderr.log"],
    ["bwa/log", "_bwa.log"],
  ]

  remove_bwa_refine_files = [
    ["bwa_refine/result", ".rmdup.*"],
    ["bwa_refine/log", "_rf.log"],
    ["bwa_refine_bam_validation/result", ".txt"],
    ["bwa_refine_bam_validation/log", "_o2o.log"],
  ]



  #clean_by_cutadapt(logger, project_dir=args.input, remove_patterns=remove_cutadapt_files+remove_posttrim_fastqc_files+remove_posttrim_fastq_len_files+remove_bwa_files+remove_bwa_refine_files,remove=not args.list)
  #clean_by_posttrim_fastqc(logger, project_dir=args.input, remove_patterns=remove_posttrim_fastqc_files,remove=not args.list)
  clean_bwa(logger, project_dir=args.input, remove_patterns=remove_bwa_files+remove_bwa_refine_files, remove=not args.list)
  #clean_bwa_refine(logger, project_dir=args.input, remove_patterns=remove_bwa_refine_files, remove=not args.list)

if __name__ == "__main__":
    main()
