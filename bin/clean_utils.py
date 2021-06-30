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


def find_error_samples_by_keyword_in_log(logger, log_folder, keyword, suffix):
  files=[y for y in glob(os.path.join(log_folder, "*" + suffix))]
  result = []
  for cfile in files:
    with open(cfile, "rt") as fin:
      for line in fin:
        if keyword in line:
          result.append(os.path.basename(cfile).replace(suffix,""))
          break
  return(result)

def find_error_samples_by_fastqc(logger, post_trim_fastqc_log_folder):
  return(find_error_samples_by_keyword_in_log(logger, post_trim_fastqc_log_folder, 'SequenceFormatException', "_fq.log"))

def find_error_samples_by_mutect(logger, mutect_log_folder):
  return(find_error_samples_by_keyword_in_log(logger, mutect_log_folder, 'A USER ERROR has occurred', "_mt.log"))

def find_error_samples_by_mutect2(logger, mutect2_log_folder):
  return(find_error_samples_by_keyword_in_log(logger, mutect2_log_folder, 'Exception', "_mt2.log"))

def find_error_samples_by_fastq_len(logger, fastq_len_result_folder):
  files=[y for y in glob(os.path.join(fastq_len_result_folder, "*.len.error"))]
  result = []
  for cfile in files:
    result.append(os.path.basename(cfile).replace(".len.error",""))
  return(result)

def find_error(line, errors):
  for error in errors:
    if error in line:
      return(True)
  return(False)

def find_error_samples_by_bwa(logger, bwa_dir):
  errors = ["[mem_sam_pe]", 'error', '[W::bseq_read]']

  files=[y for y in glob(os.path.join(bwa_dir, "result", "*.unsorted.bam.failed"))]
  result = set()
  for cfile in files:
    result.add(os.path.basename(cfile).replace(".unsorted.bam.failed",""))

  files=[y for y in glob(os.path.join(bwa_dir, "log", "*.log"))]
  for cfile in files:
    sample_name = os.path.basename(cfile).replace("_bwa.log","")
    if sample_name in result:
      continue

    with open(cfile, "rt") as fin:
      bFindSucceed = False
      for line in fin:
        if find_error(line, errors):
          result.add(sample_name)
          break
        
        if "[main] Real time" in line:
          bFindSucceed = True
      
      if not bFindSucceed:
        if not sample_name in result:
          result.add(sample_name)

  result = sorted(result)
  return(result)

def find_error_samples_by_bwa_refine(logger, dir):
  files=[y for y in glob(os.path.join(dir, "log", "*.log"))]
  result = []
  for cfile in files:
    with open(cfile, "rt") as fin:
      for line in fin:
        if ('CANCELLED' in line) or ("[W::bseq_read]" in line):
          result.append(os.path.basename(cfile).replace("_rf.log",""))
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
