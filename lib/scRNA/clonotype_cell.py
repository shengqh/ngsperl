import argparse
import logging
import os
import sys
import re
import json

def initialize_logger(logfile, args):
  logger = logging.getLogger('clonotype_cell')
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

def convert(json_file, output_file, logger):
  logger.info("reading %s" % json_file)
  with open(output_file, "wt") as fout:
    fout.write("barcode\tclonetype\n")
    barcodes = set()
    with open(json_file, "rt") as fin:
      data = json.load(fin)
      for record in data:
        if record["is_cell"]:
          if not record['barcode'] in barcodes:
            fout.write("%s\t%s\n" % (record['barcode'], record['info']['raw_clonotype_id']))
            barcodes.add(record['barcode'])
  logger.info("written to %s" % output_file)

def main():
  parser = argparse.ArgumentParser(description="Extract cell of clonotype",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = True
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help="Input clone type json file", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output clone type cell file")
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  
  if DEBUG:
    args.input="/data/h_vangard_1/alexander_gelbard_data/AG_3669_10X_cellranger4/VDJ/3669-AG-6/all_contig_annotations.json"
    args.output="/data/h_vangard_1/alexander_gelbard_data/AG_3669_10X_cellranger4/VDJ/3669-AG-6/all_contig_annotations.json.txt"

  check_file(args.input, parser)

  logger = initialize_logger(args.output + ".log", args)
  convert(args.input, args.output, logger)
  
if __name__ == "__main__":
    main()
