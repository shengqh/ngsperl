import os
import os.path
import errno
import sys
import shutil
import logging
import subprocess
import argparse

def initialize_logger(logfile, name):
  logger = logging.getLogger(name)
  loglevel = logging.INFO
  logger.setLevel(loglevel)

  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')    
 
  # create console handler and set level to info
  handler = logging.StreamHandler()
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  # create error file handler and set level to error
  handler = logging.FileHandler(logfile, "w")
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  return(logger)

def summary(logger, input_file, sample_name, thread, output_file):
  header_file = output_file + ".header"
  subprocess.call(f"samtools view -H {input_file} > {header_file}", shell=True)

  headers = []
  with open(header_file, "rt") as fin:
    bfirst = True
    for line in fin:
      if line.startswith("@RG") and ("SM:" in line):
        if bfirst:
          line = f"@RG\tID:1\tSM:{sample_name}\tPU:{sample_name}\tLB:{sample_name}\tPL:ILLUMINA"
          headers.append(line)
          bfirst = False
        continue
      elif line.startswith("@PG") and ("ID:samtools" in line) and ("PP:samtools" in line):
        continue
      else:
        headers.append(line.rstrip())
  
  with open(header_file, "wt") as fout:
    for line in headers:
      fout.write(line + "\n")
  
  cmd = f"samtools reheader -P {header_file} {input_file} > {output_file}"
  print(cmd)
  subprocess.call(cmd, shell=True)
  
  cmd = f"sambamba index -t {thread} {output_file}"
  print(cmd)
  subprocess.call(cmd, shell=True)

  logger.info("done")

if __name__ == '__main__':
  DEBUG=False
  NotDEBUG=not DEBUG

  parser = argparse.ArgumentParser(description="Replace RG in bam file",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bam file', required=NotDEBUG)
  parser.add_argument('-n', '--name', action='store', nargs='?', help="Input sample name", required=NotDEBUG)
  parser.add_argument('-t', '--thread', action='store', nargs='?', default="8", help="Input thread number", required=NotDEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output bam file", required=NotDEBUG)

  args = parser.parse_args()
  if DEBUG:
    args.input = "/scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/bam_download/reheader/1.slim.bam"
    args.name = "1"
    args.output = "/scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/bam_download/reheader/1.slim.rg.bam"
  
  logger = initialize_logger(args.output + ".log", "replace_rg")
  
  summary(logger, args.input, args.name, args.thread, args.output)
