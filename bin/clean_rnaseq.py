#!/usr/bin/env python3

import os
import os.path
import argparse
import sys
import logging
import re
 
from glob import glob

def run_cmd(cmd):
  print(cmd)
  os.system(cmd)

def remove_files(logger, project_dir, samples, remove_patterns):
  for sample in samples:
    for rfile in remove_patterns:
      remove_res_file = os.path.join(project_dir, rfile[0], rfile[1], "%s%s" % (sample, rfile[2] ))
      run_cmd('rm -rf %s' % remove_res_file)

def find_error_samples_by_fastqc(logger, post_trim_fastqc_log_folder):
  files=[y for y in glob(os.path.join(post_trim_fastqc_log_folder, "*.log"))]
  result = []
  for cfile in files:
    with open(cfile, "rt") as fin:
      for line in fin:
        if 'SequenceFormatException' in line:
          result.append(os.path.basename(cfile).replace("_fq.log",""))
          break
  return(result)

def clean_by_posttrim_fastqc(logger, project_dir, remove_patterns):
  post_trim_fastqc_log_folder = os.path.join(project_dir, "fastqc_post_trim", "log")
  samples = find_error_samples_by_fastqc(logger, post_trim_fastqc_log_folder)
  remove_files(logger, project_dir, samples, remove_patterns)

def find_error_samples_by_star(logger, star_dir):
  remove_patterns=["__STARtmp", "_Aligned.out.bam"]
  files=[y for y in glob(os.path.join(star_dir, "result", "*"))]
  check_suffix=["_SJ.out.tab", ".bamstat", ".count", ".splicing.bed"]
  result =set()
  for cfile in files:
    if cfile.endswith(".tmp"):
      run_cmd('rm -rf %s' % cfile)
      next

    for suffix in check_suffix:
      if cfile.endswith(suffix):
        sample = os.path.basename(cfile).replace(suffix, "")
        bamfile = os.path.join(star_dir, "result", "%s_Aligned.sortedByCoord.out.bam" % sample)
        if not os.path.exists(bamfile):
          result.add(sample)
      
    for rpattern in remove_patterns:
      if cfile.endswith(rpattern):
        sample = os.path.basename(cfile).replace(rpattern, "")
        result.add(sample)
    
  return(result)

def clean_star(logger, project_dir, remove_patterns):
  star_folder =  os.path.join(project_dir, "star_featurecount")
  samples = find_error_samples_by_star(logger, star_folder)
  #print(samples)
  remove_files(logger, project_dir, samples, remove_patterns)

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
    ["cutadapt", "result", "_clipped*"],
    ["cutadapt", "result", ".version"],
    ["cutadapt", "log", "_cut.log"],
    ["fastq_len", "result", ".len"],
    ["fastq_len", "log", "_flen.log"],
  ]

  remove_posttrim_fastqc_files = [
    ["fastqc_post_trim", "result", "/*"],
    ["fastqc_post_trim", "log", "_fq.log"],
  ]

  remove_star_files = [
    ["star_featurecount", "result", ".count*"],
    ["star_featurecount", "result", ".bamstat"],
    ["star_featurecount", "result", "__STARtmp"],
    ["star_featurecount", "result", "_Aligned*"],
    ["star_featurecount", "result", "_Log*"],
    ["star_featurecount", "result", "_SJ.out.tab"],
    ["star_featurecount", "result", ".splicing.bed"],
    ["star_featurecount", "log", "_sf.log"],
  ]

  clean_by_posttrim_fastqc(logger, project_dir=args.input, remove_patterns=remove_cutadapt_files+remove_posttrim_fastqc_files+remove_star_files)
  clean_star(logger, project_dir=args.input, remove_patterns=remove_star_files)
  
if __name__ == "__main__":
    main()
