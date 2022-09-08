#!/usr/bin/env python3

import os
import os.path
import argparse
import sys
import logging
import re
import subprocess
 
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


def get_sample_map(files, file_pattern, name_pattern, use_dir_name=0, recursive_dir=False, fill_zero=False, add_prefix_P=False, add_prefix_X=False):
  pattern = re.compile(file_pattern)
  files = sorted([f for f in files if pattern.search(os.path.basename(f))])

  #print(files)

  result={}
  if use_dir_name > 0:
    for f in files:
      sample_name = get_dir_name(f, use_dir_name)
      sample_name = get_sample_name(sample_name, name_pattern)
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
  return(result)

def find_files(logger, source_dir, output_file, file1_pattern, file2_pattern, name_pattern, use_dir_name=0, recursive_dir=False, fill_zero=False, add_prefix_P=False, add_prefix_X=False):
  if recursive_dir:
    output = subprocess.getoutput(f"gsutil ls -r {source_dir}/**")
    print(output)
  else:
    output = subprocess.getoutput(f"gsutil ls {source_dir}")
  files = output.split('\n')

  read1map = get_sample_map(files, file1_pattern, name_pattern, use_dir_name, recursive_dir, fill_zero, add_prefix_P, add_prefix_X)
  read2map = get_sample_map(files, file2_pattern, name_pattern, use_dir_name, recursive_dir, fill_zero, add_prefix_P, add_prefix_X)

  prefix = "P" if add_prefix_P else ""
  prefix = "X" if add_prefix_X else prefix

  with open(output_file, "wt") as fout:
    fout.write("entity:fastq_id\tfastq_1\tfastq_2\n")
    for sample_name in sorted(read1map.keys()):
      read1_files = read1map[sample_name]
      read2_files = read2map[sample_name]
      read1_files_str = '", "'.join(read1_files)
      read2_files_str = '", "'.join(read2_files)
      o_sample_name = prefix + sample_name.replace("-", "_")
      fout.write(f"\"{o_sample_name}\"\t[\"{read1_files_str}\"]\t[\"{read2_files_str}\"]\n")
        
def main():
  parser = argparse.ArgumentParser(description="Get gcp file definition",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input source folder', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help='Output terra data file', required=NOT_DEBUG)
  parser.add_argument('--read1_pattern', action='store', nargs='?', help='Input read1 pattern', required=NOT_DEBUG)
  parser.add_argument('--read2_pattern', action='store', nargs='?', help='Input read2 pattern', required=NOT_DEBUG)
  parser.add_argument('-n', '--name_pattern', action='store', nargs='?', help='Input name pattern')
  parser.add_argument('-a', '--add_zero', action='store_true', help='Filling digits with zero')
  parser.add_argument('-d', '--use_dir_name', action='store', type=int, default=0, help='Use X level dir name as sample name, 0 means use file name')
  parser.add_argument('-r', '--recursive_dir', action='store_true', help='Find file in recursive dir')
  parser.add_argument('-p', '--add_prefix_P', action='store_true', help='Add P as prefix of sample name')
  parser.add_argument('-x', '--add_prefix_X', action='store_true', help='Add X as prefix of sample name')
  
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

  if not args.input.startswith("gs://"):
    args.input= os.path.abspath(args.input)
  print(args)
  
  logger = logging.getLogger('file_def')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  find_files(logger, source_dir=args.input, output_file=args.output, file1_pattern=args.read1_pattern, file2_pattern=args.read2_pattern, name_pattern=args.name_pattern,
             use_dir_name=args.use_dir_name, recursive_dir=args.recursive_dir, fill_zero=args.add_zero, 
             add_prefix_P=args.add_prefix_P, add_prefix_X=args.add_prefix_X)
  
if __name__ == "__main__":
    main()
