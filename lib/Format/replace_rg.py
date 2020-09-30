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

def reheader(logger, input_file, sample_name, sorted, replace_PL, thread, output_file):
  header_file = output_file + ".header"
  subprocess.call(f"samtools view -H {input_file} > {header_file}", shell=True)

  ids = set()
  headers = []
  bFindRG = False
  with open(header_file, "rt") as fin:
    for line in fin:
      if line.startswith("@HD") and ("SO:" in line):
        line = line.rstrip()
        if sorted != None:
          line = line.split("SO:")[0] + "SO:" + sorted
        headers.append(line)
        continue

      # if line.startswith("@SQ") and ("SN:pBACe3.6" in line):
      #   line = "@SQ\tSN:pBACe3.6\tLN:11612\tAS:combined_aav_genome.ndx"
      #   headers.append(line)
      #have to keep old ID, PU and LB, otherwise the field in read "RG:Z" (id) will be wrong
      if line.startswith("@RG") and ("SM:" in line):
        bFindRG = True
        parts = line.rstrip().split('\t')
        for idx in range(0, len(parts)):
          part = parts[idx]
          if part.startswith("SM"):
            parts[idx] = f"SM:{sample_name}"
          if part.startswith("PL") and replace_PL:
            parts[idx] = f"PL:ILLUMINA"
        if not "PU:" in line:
          parts.append(f"PU:{sample_name}")
        if not "LB:" in line:
          parts.append(f"LB:{sample_name}")
        if not "PL:" in line:
          parts.append("PL:ILLUMINA")
        line = "\t".join(parts)
        headers.append(line)
        continue
      elif line.startswith("@PG") and ("ID:samtools" in line) and (("PP:samtools" in line) or ("PN:samtools" in line)):
        continue
      else:
        headers.append(line.rstrip())
  
  if not bFindRG:
    oldHeaders = headers
    headers = []
    bSQ = False
    bAdded = False
    for line in oldHeaders:
      if line.startswith("@SQ"):
        bSQ = True
        headers.append(line)
      elif not bSQ:
        headers.append(line)
      elif not bAdded:
        headers.append(f"@RG\tID:1\tSM:{sample_name}\tPL:ILLUMINA");
        bAdded = True
        headers.append(line)
      else:
        headers.append(line)

  with open(header_file, "wt") as fout:
    for line in headers:
      fout.write(line + "\n")
  
  tmpfile = output_file + ".tmp.bam"
  cmd = f"samtools reheader -P {header_file} {input_file} > {tmpfile}"
  print(cmd)
  subprocess.call(cmd, shell=True)
  
  cmd = f"sambamba index -t {thread} {tmpfile}"
  print(cmd)
  subprocess.call(cmd, shell=True)

  if os.path.isfile(output_file):
    os.remove(output_file)
  os.rename(tmpfile, output_file)

  if os.path.isfile(output_file + ".bai"):
    os.remove(output_file + ".bai")
  os.rename(tmpfile + ".bai", output_file + ".bai")

  logger.info("done")

if __name__ == '__main__':
  DEBUG=False
  NotDEBUG=not DEBUG

  parser = argparse.ArgumentParser(description="Replace RG in bam file",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bam file', required=NotDEBUG)
  parser.add_argument('-n', '--name', action='store', nargs='?', help="Input sample name", required=NotDEBUG)
  parser.add_argument('-t', '--thread', action='store', nargs='?', default="8", help="Input thread number", required=NotDEBUG)
  parser.add_argument('-s', '--sorted', action='store', nargs='?', choices=["coordinate", "queryname"], help="Replace with sorted (coordinate or queryname, default as original)")
  parser.add_argument('--replace_PL', action='store_true', help="Replace PL with ILLUMINA")
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output bam file", required=NotDEBUG)

  args = parser.parse_args()
  if DEBUG:
    args.input = "/scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/bam_download/reheader/1.slim.bam"
    args.name = "1"
    args.output = "/scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/bam_download/reheader/1.slim.rg.bam"
  
  logger = initialize_logger(args.output + ".log", "replace_rg")
  
  reheader(logger, args.input, args.name, args.sorted, args.replace_PL, args.thread, args.output)
