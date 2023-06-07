import argparse
import logging
import os
import sys
import gzip

def initialize_logger(logfile, args):
  logger = logging.getLogger('fragment_cells')
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

def get_cells(input_file, output_file, logger):
  logger.info(f"processing {input_file} ...")

  cell_hash = {}
  count = 0
  with gzip.open(input_file, "rt") as fin:
    for line in fin:
      if line.startswith('#'):
        continue
      parts = line.split('\t')
      cell = parts[3]
      cell_hash[cell] = 1
      count += 1
      if count % 1000000 == 0:
        logger.info(f"{count} : {len(cell_hash)} cells")

  all_cells = sorted(list(cell_hash.keys()))
  with open(output_file, "wt") as fout:
    for cell in all_cells:
      fout.write(f"{cell}\n")

  logger.info("done")

def main():
  parser = argparse.ArgumentParser(description="get fragment cells",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help="Input fragment file", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output fragment cell file")
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  
  if DEBUG:
    args.input="/data/jbrown_lab/2023/20230503_encode_scRNA_snATAC/ENCDO068KYD/encode_scatac_dcc_2/results/ENCSR650BBI-1/fragments/fragments.tsv.gz"
    args.output="/data/jbrown_lab/2023/20230503_encode_scRNA_snATAC/ENCDO068KYD/encode_scatac_dcc_2/results/ENCSR650BBI-1/fragments/fragments.cells.txt"

  logger = initialize_logger(args.output + ".log", args)
  get_cells(args.input, args.output, logger)
  
if __name__ == "__main__":
    main()
