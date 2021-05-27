#!/usr/bin/env python3

import os
import os.path
import argparse
import sys
import logging
import re
 
from glob import glob

def get_dir_name(file_name, use_dir_name):
  result = file_name
  while use_dir_name > 0:
    result=os.path.dirname(result)
    use_dir_name -= 1
  result=os.path.basename(result)
  return(result)

def get_sample_name(f, name_pattern):
  result = os.path.basename(f)
  if name_pattern != None:
    m = re.search(name_pattern, result)
    if m != None:
      result=m.group()
  return(result)

def find_files(logger, source_dir, file_pattern, name_pattern, use_dir_name=0, recursive_dir=False):
  if recursive_dir:
    files=[y for x in os.walk(source_dir) for y in glob(os.path.join(x[0], file_pattern))]
  else:
    files=[y for y in glob(os.path.join(source_dir, file_pattern))]

  #print(files)

  result={}
  if use_dir_name > 0:
    for f in files:
      sample_name = get_dir_name(f, use_dir_name)
      if sample_name not in result:
        result[sample_name] = [f]
      else:
        result[sample_name].append(f)
  else:
    for f in files:
      sample_name = get_sample_name(f, name_pattern)
      if sample_name not in result:
        result[sample_name] = [f]
      else:
        result[sample_name].append(f)
  
  print("  files => {")
  for sample_name in sorted(result.keys()):
    sample_files = result[sample_name]
    sample_files = sorted(sample_files)
    print("    '%s' => [ '%s' ], " % (sample_name, "', '".join(sample_files) ))
  print("  }, ")

  return(result)
        
def main():
  parser = argparse.ArgumentParser(description="Get file definition",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input source folder', required=NOT_DEBUG)
  parser.add_argument('-f', '--file_pattern', action='store', nargs='?', help='Input file pattern', required=NOT_DEBUG)
  parser.add_argument('-n', '--name_pattern', action='store', nargs='?', help='Input name pattern')
  parser.add_argument('-d', '--use_dir_name', action='store', type=int, default=0, help='Use X level dir name as sample name, 0 means use file name')
  parser.add_argument('-r', '--recursive_dir', action='store_true', help='Find file in recursive dir')
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    args.input="/scratch/cqs/shengq2/alexander_gelbard_projects/20210219_scRNA_JLin_human/VDJ_chain/result"
    #args.input="/scratch/cqs/shengq2/alexander_gelbard_projects/20210219_scRNA_JLin_human/VDJ_chain/result/nMCu81420/outs"
    args.file_pattern="all_contig_annotations.json"
    args.use_dir_name = 2
    args.recursive_dir = True

  args.input= os.path.abspath(args.input)
  print(args)
  
  logger = logging.getLogger('file_def')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  find_files(logger, source_dir=args.input, file_pattern=args.file_pattern, name_pattern=args.name_pattern,
             use_dir_name=args.use_dir_name, recursive_dir=args.recursive_dir)
  
if __name__ == "__main__":
    main()
