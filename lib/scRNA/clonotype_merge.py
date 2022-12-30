import argparse
import logging
import os
import sys
import re
import json
from collections import OrderedDict

def initialize_logger(logfile, args):
  logger = logging.getLogger('clonotype_merge')
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

def getValidFilename(s):
  s = str(s).strip().replace(' ', '_')
  return re.sub(r'(?u)[^-\w.]', '', s)
  
def runCommand(command, logger):
  logger.info("run : " + command )
  os.system(command)
  
def check_file(filename, parser):
  if not os. path. isfile(filename):
    print("error: file not exists: " + filename)
    parser.print_help()
    sys.exit(1)

def read_file_map(fileName):
  result = OrderedDict()
  with open(fileName) as fh:
    for line in fh:
      filepath, name = line.strip().split('\t', 1)
      result[name] = filepath.strip()
  return(result)

def merge(json_file_list, output_file, logger):
  json_files = read_file_map(json_file_list)
  finaldata = []
  cdr3map = {}
  for json_file_name in json_files:
    json_file = json_files[json_file_name]
    logger.info("reading %s" % json_file)
    with open(json_file, "rt") as fin:
      data = json.load(fin)
      for record in data:
        record['barcode'] = json_file_name + "_" + record['barcode']
        if 'info' in record:
          for key in record['info']:
            v = record['info'][key]
            if v != None:
              record['info'][key] = json_file_name + "_" + record['info'][key]
        finaldata.append(record)

        cdr3=record['cdr3']
        if cdr3 is None:
          continue

        for ann in record['annotations']:
          chain = ann['feature']['chain']
          cdr3map[cdr3] = chain
  
  logger.info("writing %s" % output_file)
  with open(output_file, "wt") as fout:
    json.dump(finaldata, fout, indent=4)
  
  with open(output_file + ".cdr3", "wt") as fout:
    fout.write("cdr3\tchain\n")
    for cdr3 in sorted(cdr3map.keys()):
      fout.write(f"{cdr3}\t{cdr3map[cdr3]}\n")

  logger.info("done")

def main():
  parser = argparse.ArgumentParser(description="merge clonotype data",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help="Input clone type json file list", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output folder")
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  
  if DEBUG:
    args.input="/scratch/cqs/alexander_gelbard_projects/20201117_scRNA_3669_enclone/clonotype_merge/result/AG3669__fileList1.list"
    args.output="/scratch/cqs/alexander_gelbard_projects/20201117_scRNA_3669_enclone/clonotype_merge/result/all_contig_annotations.json"

  check_file(args.input, parser)

  logger = initialize_logger(args.output + ".log", args)
  merge(args.input, args.output, logger)
  
if __name__ == "__main__":
    main()
