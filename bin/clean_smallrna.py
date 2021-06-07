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
  cutadapt_result_folder = os.path.join(project_dir, "intermediate_data", "cutadapt", "result")
  samples = do_clean_by_cutadapt(logger, project_dir, cutadapt_result_folder, remove_patterns, remove=remove)
  return(samples)

def clean_fastq_len(logger, project_dir, remove_patterns, remove=True):
  fastq_len_result_folder = os.path.join(project_dir, "preprocessing", "fastq_len", "result")
  samples = find_error_samples_by_fastq_len(logger, fastq_len_result_folder)
  remove_files(logger, project_dir, samples, remove_patterns, remove)
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
    ["intermediate_data/cutadapt/result", "_clipped*"],
    ["intermediate_data/cutadapt/result", ".version"],
    ["intermediate_data/cutadapt/log", "_cut.log"],
  ]

  remove_posttrim_fastqc_files = [
    ["intermediate_data/fastqc_post_trim/result", "/*"],
    ["intermediate_data/fastqc_post_trim/log", "_fq.log"],
  ]

  remove_posttrim_fastq_len_files = [
    ["preprocessing/fastq_len/result", ".len"],
    ["preprocessing/fastq_len/result", ".len.error"],
  ]

  remove_posttrim_files = [
    ["intermediate_data/bowtie1*/result", ".bam*"],
    ["intermediate_data/bowtie1_genome_1mm_NTA_pmnames/result", ".pmnames"],
    ["intermediate_data/bowtie1*count/result", "/*"],
    ["preprocessing/identical/result", "_clipped*"],
    ["intermediate_data/identical_NTA/result","_clipped*"],
    ["intermediate_data/identical_check_cca/result","_CCA.tsv"]
  ]

  remove=True
  clean_by_cutadapt(logger, project_dir=args.input, remove_patterns=remove_cutadapt_files+remove_posttrim_fastqc_files+remove_posttrim_fastq_len_files+remove_posttrim_files, remove=remove)
  clean_fastq_len(logger, project_dir=args.input, remove_patterns=remove_cutadapt_files+remove_posttrim_fastqc_files+remove_posttrim_fastq_len_files+remove_posttrim_files, remove=remove)

if __name__ == "__main__":
    main()
