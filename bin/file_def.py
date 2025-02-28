#!/usr/bin/env python3

import os
import os.path
import argparse
import sys
import logging
import re
import csv
 
from glob import glob

def get_valid_sample_name(old):
  return re.sub('[ \-*+]+', '_', old)

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
      result="".join(m.groups())
  return(result)

def read_meta(meta_file, meta_file_index, meta_name_index, meta_file_delimiter=','):
  result = {}
  with open(meta_file, newline='') as f:
    reader = csv.reader(f, delimiter=meta_file_delimiter)
    for row in reader:  
      result[row[meta_file_index]] = row[meta_name_index]
  return(result)

def replace_name(old, func):
  result = {}
  for sample_name in sorted(old.keys()):
    new_name = func(sample_name)
    result[new_name] = old[sample_name]
  return(result)

def find_files(logger, source_dir, file_pattern, name_pattern, use_dir_name=0, recursive_dir=False, fill_zero=False, add_prefix_P=False, add_prefix_X=False, 
  meta_file=None, meta_file_index=-1, meta_name_index=-1,meta_file_delimiter=','):

  if meta_file != None:
    meta_dic = read_meta(meta_file, meta_file_index, meta_name_index)
    print("  meta data:")
    for key in sorted(meta_dic.keys()):
      print(f"    {key} : {meta_dic[key]}")
    print("")
  else:
    meta_dic = None

  if recursive_dir:
    files=[y for x in os.walk(source_dir) for y in glob(os.path.join(x[0], file_pattern))]
  else:
    files=[y for y in glob(os.path.join(source_dir, file_pattern))]

  #print(files)

  result={}
  if use_dir_name > 0:
    for f in files:
      f=f.replace("/panfs/accrepfs.vampire", "")
      sample_name = get_dir_name(f, use_dir_name)
      sample_name = get_sample_name(sample_name, name_pattern)
      if sample_name not in result:
        result[sample_name] = [f]
      else:
        result[sample_name].append(f)
  else:
    for f in files:
      f=f.replace("/panfs/accrepfs.vampire", "")
      sample_name = get_sample_name(f, name_pattern)
      if sample_name not in result:
        result[sample_name] = [f]
      else:
        result[sample_name].append(f)
  
  if meta_file != None:
    result = replace_name(result, lambda x: meta_dic[x] if x in meta_dic else x)

  if fill_zero:
    numbers = {}
    max_len = 0
    for sample_name in sorted(result.keys()):
      digits = re.search('(\d+)$', sample_name)
      if digits:
        numbers[sample_name] = digits.group(1)
        max_len = max(max_len, len(numbers[sample_name]))
        #print("%s => %s" % (sample_name, numbers[sample_name]))
    
    new_result = {}
    for sample_name in result.keys():
      if sample_name not in numbers:
        new_result[sample_name] = result[sample_name]
        continue
      
      number = numbers[sample_name]
      if len(number) == max_len:
        new_result[sample_name] = result[sample_name]
        continue

      new_number = number.zfill(max_len)
      new_name = re.sub('\d+$', new_number, sample_name)
      #print("%s : %s => %s => %s" % (sample_name, number, new_number, new_name))
      new_result[new_name] = result[sample_name]

    result = new_result

  prefix = "P" if add_prefix_P else ""
  prefix = "X" if add_prefix_X else prefix

  result = replace_name(result, lambda x: prefix+get_valid_sample_name(x))

  print("  files => {")
  for sample_name in sorted(result.keys()):
    sample_files = sorted(result[sample_name])
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
  parser.add_argument('-a', '--add_zero', action='store_true', help='Filling digits with zero')
  parser.add_argument('-d', '--use_dir_name', action='store', type=int, default=0, help='Use X level dir name as sample name, 0 means use file name')
  parser.add_argument('-r', '--recursive_dir', action='store_true', help='Find file in recursive dir')
  parser.add_argument('-p', '--add_prefix_P', action='store_true', help='Add P as prefix of sample name')
  parser.add_argument('-x', '--add_prefix_X', action='store_true', help='Add X as prefix of sample name')
  parser.add_argument('--meta_file', action='store', nargs='?', help='Input meta file')
  parser.add_argument('--meta_file_index', action='store', type=int, help='Input index corresponding to real file name')
  parser.add_argument('--meta_name_index', action='store', type=int, help='Input index corresponding to expected file name')
  parser.add_argument('--meta_file_delimiter', action='store', nargs='?', help='Input meta file delimiter')
  
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
             use_dir_name=args.use_dir_name, recursive_dir=args.recursive_dir, fill_zero=args.add_zero, 
             add_prefix_P=args.add_prefix_P, add_prefix_X=args.add_prefix_X,
             meta_file=args.meta_file, meta_file_index=args.meta_file_index, meta_name_index=args.meta_name_index, meta_file_delimiter=args.meta_file_delimiter)
  
if __name__ == "__main__":
    main()
