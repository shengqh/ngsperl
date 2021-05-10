#!/usr/bin/env python3

import os
import os.path
import logging
 
from glob import glob

def run_cmd(cmd, remove=True):
  print(cmd)
  if remove:
    os.system(cmd)

def remove_files(logger, project_dir, samples, remove_patterns, remove=True):
  for sample in samples:
    for rfile in remove_patterns:
      remove_res_file = project_dir + "/" + rfile[0] + "/" + sample + rfile[1]
      run_cmd('rm -rf %s' % remove_res_file, remove)

def find_error_samples_by_cutadapt(logger, cutadapt_result_folder):
  files=[y for y in glob(os.path.join(cutadapt_result_folder, "*.fastq.gz"))]
  result = []
  for cfile in files:
    sfile=cfile.split("_clipped")[0]
    versionfile=sfile + ".version"
    if not os.path.exists(versionfile):
      result.append(os.path.basename(sfile))
  return(result)

def do_clean_by_cutadapt(logger, project_dir, cutadapt_result_folder, remove_patterns, remove=True):
  samples = find_error_samples_by_cutadapt(logger, cutadapt_result_folder)
  remove_files(logger, project_dir, samples, remove_patterns, remove=remove)
  return(samples)

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

def find_error_samples_by_fastq_len(logger, fastq_len_result_folder):
  files=[y for y in glob(os.path.join(fastq_len_result_folder, "*.len.error"))]
  result = []
  for cfile in files:
    result.append(os.path.basename(cfile).replace(".len.error",""))
  return(result)

def find_error_samples_by_bwa(logger, bwa_dir):
  files=[y for y in glob(os.path.join(bwa_dir, "log", "*.log"))]
  result = []
  for cfile in files:
    with open(cfile, "rt") as fin:
      for line in fin:
        if 'error' in line:
          result.append(os.path.basename(cfile).replace("_bwa.log",""))
          break
  return(result)

def find_error_samples_by_star(logger, star_dir):
  remove_patterns=[".0000.bam", "__STARgenome", "__STARtmp", "_Aligned.out.bam"]
  files=[y for y in glob(os.path.join(star_dir, "result", "*"))]
  check_suffix=["_SJ.out.tab", ".bamstat", ".count", ".splicing.bed"]
  result =set()
  for cfile in files:
    if cfile.endswith(".tmp"):
      run_cmd('rm -rf %s' % cfile)
      next

    for suffix in check_suffix:
      if cfile.endswith(suffix) and not cfile.endswith(".chromosome.count"):
        sample = os.path.basename(cfile).replace(suffix, "")
        bamfile = os.path.join(star_dir, "result", "%s_Aligned.sortedByCoord.out.bam" % sample)
        if not os.path.exists(bamfile):
          result.add(sample)
      
    for rpattern in remove_patterns:
      if cfile.endswith(rpattern):
        sample = os.path.basename(cfile).replace(rpattern, "")
        result.add(sample)
  result=sorted(result) 
  return(result)
