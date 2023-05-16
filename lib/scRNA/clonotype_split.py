import argparse
import logging
import os
import os.path
import sys
import re
import json
import pandas as pd
from collections import OrderedDict

def initialize_logger(logfile, args):
  logger = logging.getLogger('clonotype_split')
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
  
def check_file(filename, parser):
  if not os. path. isfile(filename):
    print("error: file not exists: " + filename)
    parser.print_help()
    sys.exit(1)

def split(json_file, cell_hashtag_file, hashtag_sample_file, output_folder, logger):
  logger.info("reading %s" % cell_hashtag_file)
  cells=pd.read_csv(cell_hashtag_file)
  cells=cells.loc[cells['HTO.global'] == 'Singlet']
  barcode_dict = dict(zip(cells.iloc[:, 0], cells.HTO))
  #print(barcode_dict)

  logger.info("reading %s" % hashtag_sample_file)
  samples=pd.read_table(hashtag_sample_file, header=None)
  samples_dict = dict(zip(samples.iloc[:, 1], samples.iloc[:, 0]))
  print(samples_dict)

  logger.info("reading %s" % json_file)
  json_data = []
  with open(json_file, "rt") as fin:
    data = json.load(fin)
    for record in data:
      json_data.append(record)
  
  for sample_name in samples_dict.values():
    sample_folder = os.path.join(output_folder, sample_name)
    if not os.path.isdir(sample_folder):
      os.mkdir(sample_folder)
    
    sample_file = os.path.join(sample_folder, "all_contig_annotations.json")
    sample_data = []
    for record in json_data:
      if record['barcode'] in barcode_dict:
        if barcode_dict[record['barcode']] in samples_dict:
          if samples_dict[barcode_dict[record['barcode']]] == sample_name:
            sample_data.append(record)
    logger.info("writing %s" % sample_file)
    with open(sample_file, "wt") as fout:
      json.dump(sample_data, fout, indent=4)

  logger.info("done")

def main():
  parser = argparse.ArgumentParser(description="merge clonotype data",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help="Input clone type json file", required=NOT_DEBUG)
  parser.add_argument('-c', '--cell_hashtag', action='store', nargs='?', help="Input cell hashtag file", required=NOT_DEBUG)
  parser.add_argument('-s', '--hashtag_sample', action='store', nargs='?', help="Input hashtag sample file", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output folder")
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  
  if DEBUG:
    args.input="/data/cqs/annet_kirabo_data/20210809_6383_AP_ECCITE/TCR/6383-AP-1/all_contig_annotations.json"
    args.cell_hashtag="/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/hto_samples_cutoff_souporcell_integration/result/AP_1.HTO.csv"
    args.hashtag_sample="/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/clonotype_01_split/result/fileList_3_AP_1.txt"
    args.output="/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/clonotype_01_split/result/"

  check_file(args.input, parser)
  check_file(args.cell_hashtag, parser)
  check_file(args.hashtag_sample, parser)

  logger = initialize_logger(os.path.join(args.output, "clonotype_split.log"), args)
  split(args.input, args.cell_hashtag, args.hashtag_sample, args.output, logger)
  
if __name__ == "__main__":
    main()
