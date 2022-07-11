import argparse
import logging
import os
import sys
import re
import csv
from collections import OrderedDict

def initialize_logger(logfile, args):
  logger = logging.getLogger('enclone_to_clonotype')
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

def sortByCells(val): 
    return val['n'] 

def convert(enclone_file, chain_file, output_file, logger):
  cdr3map = {}
  with open(chain_file) as fin:
    fin.readline()
    for line in fin:
      cdr3, chain = line.strip().split('\t', 1)
      cdr3map[cdr3] = chain

  rows = []
  with open(enclone_file, "rt") as fin:
    reader = csv.DictReader(fin)
    for row in reader:
      rows.append(row)

  for row in rows:
    row['n'] = int(row['n'])

  #rows.sort(key=sortByCells, reverse=True)

  clone_sample_cell_map = {}
  samples = []
  clono_index = 0
  for row in rows:
    clono_index += 1
    sample_map = {}
    clone_sample_cell_map[clono_index] = sample_map
    barcodes = row['barcodes'].split(',')
    if '_' in barcodes[0]:
      for barcode in barcodes:
        m = re.search('(.+)_', barcode)
        sample = m.group(1)
        samples.append(sample)
        if sample in sample_map:
          sample_map[sample] += 1
        else:
          sample_map[sample] = 1

  samples = sorted([s for s in set(samples)])
  print(samples)

  total_cells = sum(row['n'] for row in rows)
  with open(output_file, "wt") as fout:
    clono_index = 0
    fout.write("clonotype_id,frequency,proportion,cdr3s_aa,cdr3s_nt,%s,TRBV,TRBJ,cells\n" % ",".join(samples))
    for row in rows:
      clono_index += 1
      sample_map = clone_sample_cell_map[clono_index]
      aa1 = row['cdr3_aa1']
      aa2 = row['cdr3_aa2']
      dna1 = row['cdr3_dna1']
      dna2 = row['cdr3_dna2']
      trbv=row['v_name1']
      trbj=row['j_name1']
      aas = []
      dnas = []
      if aa1 != "":
        chain = cdr3map[aa1]
        aas.append(chain + ":" + aa1)
        dnas.append(chain + ":" + dna1)
      if aa2 != "":
        chain = cdr3map[aa2]
        aas.append(chain + ":" + aa2)
        dnas.append(chain + ":" + dna2)
      aas.sort()
      dnas.sort()
      barcodes = row['barcodes']
      barcodes = barcodes.replace(",",";")
      cell_counts = []
      for sample in samples:
        cell_counts.append(str(sample_map[sample]) if sample in sample_map else "0")
      fout.write("clonotype%d,%s,%s,%s,%s,%s,%s,%s,%s\n" % (
        clono_index, 
        row['n'], 
        row['n']*1.0/total_cells, 
        ";".join(aas), 
        ";".join(dnas), 
        ",".join(cell_counts),
        trbv,
        trbj,
        barcodes
        ))

  logger.info("done")

def main():
  parser = argparse.ArgumentParser(description="merge clonotype data",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help="Input enclone file", required=NOT_DEBUG)
  parser.add_argument('-c', '--chain', action='store', nargs='?', help="Input cdr3 chain file", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output clonotype file", required=NOT_DEBUG)
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  
  if DEBUG:
    args.input="/scratch/cqs/alexander_gelbard_projects/20201202_scRNA_5126/clonotype_merge_enclone/result/scRNA_5126.csv"
    args.chain="/scratch/cqs/alexander_gelbard_projects/20201202_scRNA_5126/clonotype_merge/result/all_contig_annotations.json.cdr3"
    args.output="/scratch/cqs/alexander_gelbard_projects/20201202_scRNA_5126/clonotype.csv"

  check_file(args.input, parser)
  check_file(args.chain, parser)

  logger = initialize_logger(args.output + ".log", args)
  convert(args.input, args.chain, args.output, logger)
  
if __name__ == "__main__":
    main()
