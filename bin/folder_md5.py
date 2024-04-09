#!/usr/bin/env python3

import os
import os.path
import argparse
import sys
import logging

from glob import glob
import hashlib

def folder_md5( logger, 
                source_dir, 
                file_pattern, 
                recursive_dir,
                output_file):

  if recursive_dir:
    files=[y for x in os.walk(source_dir) for y in glob(os.path.join(x[0], file_pattern))]
  else:
    files=[y for y in glob(os.path.join(source_dir, file_pattern))]

  with open(output_file, "wt") as fout:
    for f in files:
      with open(f, 'rb') as file:
        data = file.read()
        md5_hash = hashlib.md5(data).hexdigest()
        msg = f"{md5_hash}\t{os.path.basename(f)}"
        fout.write(f"{msg}\n")
        logger.info(f"{msg}")     

  logger.info(f"Output file: {output_file}")   

def main():
  parser = argparse.ArgumentParser(description="Get file definition",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input source folder', required=NOT_DEBUG)
  parser.add_argument('-f', '--file_pattern', action='store', nargs='?', help='Input file pattern', required=NOT_DEBUG)
  parser.add_argument('-r', '--recursive_dir', action='store_true', help='Find file in recursive dir')
  parser.add_argument('-o', '--output', action='store', nargs='?', help='Output file (default is the folder name as prefix)', required=False)

  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    args.input="/data/shah_lab/shengq2/20220220_obesity_smallRNA"
    args.file_pattern="*.fastq.gz"
    args.recursive_dir=True
    args.output="20220220_obesity_smallRNA.md5.txt"

  args.input= os.path.abspath(args.input)
  if(args.output == None):
    args.output = os.path.basename(args.input) + ".md5.txt"
  print(args)
  
  logger = logging.getLogger('folder_md5')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  folder_md5( logger, 
              source_dir=args.input, 
              file_pattern=args.file_pattern, 
              recursive_dir=args.recursive_dir,
              output_file=args.output)
  
if __name__ == "__main__":
    main()
