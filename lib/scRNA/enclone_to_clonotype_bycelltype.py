import argparse
import logging
import os
import sys
import re
import csv
from collections import OrderedDict
import pandas as pd

def initialize_logger():
  logger = logging.getLogger('enclone_to_clonotype')
  loglevel = logging.INFO
  logger.setLevel(loglevel)

  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')    
 
  # create console handler and set level to info
  handler = logging.StreamHandler()
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

def getValue(rows, colname):
  result = ""
  for row in rows:
    if row[colname] != "":
      result = row[colname]
      break
  return(result)

def getCdr3AA(rows, colname, cdr3map):
  result = getValue(rows, colname)
  chain = ""
  if result != "":
    chain = cdr3map[result]
    result = chain + ":" + result
  return(result, chain)

def getCdr3DNA(rows, colname, chain):
  result = getValue(rows, colname)
  if result != "":
    result = chain + ":" + result
  return(result)

def getChain(idx, rows, cdr3map):
  index = str(idx)
  aa, chain = getCdr3AA(rows, 'cdr3_aa' + index, cdr3map)
  dna = getCdr3DNA(rows, 'cdr3_dna' + index, chain)
  v_name = getValue(rows, 'v_name' + index)
  d_name = getValue(rows, 'd_name' + index)
  j_name = getValue(rows, 'j_name' + index)
  return {"aa":aa, "dna":dna, "v_name":v_name, "d_name":d_name, "j_name":j_name}

def readCsv(filename, sep=","):
  allrows = []
  with open(filename, "rt") as fin:
    reader = csv.DictReader(fin, delimiter=sep)
    for row in reader:
      allrows.append(row)
  return(allrows)

def mergeChain(clono_data, value_map):
  if value_map['aa'] != "":
    clono_data['aa'].append(value_map['aa'])
    clono_data['dna'].append(value_map['dna'])
    clono_data['v_name'].append(value_map['v_name'])
    clono_data['d_name'].append(value_map['d_name'])
    clono_data['j_name'].append(value_map['j_name'])
  return(clono_data)

def add_clono_data(cid, frequency, proportion, rows, cdr3map, samples, clono_types):
  clono_data = {"clonotype_id":cid, "frequency":frequency, "proportion":proportion, "aa":[], "dna":[], "v_name":[], "d_name":[],"j_name":[]}
  for chain in range(1,5):
    chain_map = getChain(chain, rows, cdr3map)
    clono_data = mergeChain(clono_data, chain_map)

  sample_map = {}
  all_barcodes = [row['barcode'] for row in rows]
  for barcode in all_barcodes:
    m = re.search('(.+)_', barcode)
    sample = m.group(1)
    samples.append(sample)
    if sample in sample_map:
      sample_map[sample] += 1
    else:
      sample_map[sample] = 1
  clono_data['cells'] = sample_map
  clono_data['barcodes'] = sorted(all_barcodes)
  clono_types.append(clono_data)

def write_to_file(output_file, clono_types, samples):
  with open(output_file, "wt") as fout:
    fout.write("clonotype_id,frequency,proportion,cdr3s_aa,cdr3s_nt,%s,v_names,d_names,j_names,cells\n" % ",".join(samples))
    for clono_data in clono_types:
      sample_map = clono_data['cells']
      cell_counts = []
      for sample in samples:
        cell_counts.append(str(sample_map[sample]) if sample in sample_map else "0")

      #clono_data = {"clonotype_id":cid, "frequency":frequency, "proportion":proportion, "aa":[], "dna":[], "v_name":[], "d_name":[],"j_name":[]}

      fout.write("%s,%d,%f,%s,%s,%s,%s,%s,%s,%s\n" % (
        clono_data['clonotype_id'], 
        clono_data['frequency'], 
        clono_data['proportion'], 
        ";".join(clono_data['aa']), 
        ";".join(clono_data['dna']), 
        ",".join(cell_counts),
        ";".join(clono_data['v_name']), 
        ";".join(clono_data['d_name']), 
        ";".join(clono_data['j_name']),
        ";".join(clono_data['barcodes'])
        ))

def convert(logger, enclone_file, chain_file, celltype_file, output_file_prefix):
  cdr3map_rows = readCsv(chain_file, sep="\t")
  cdr3map = {row['cdr3']:row['chain'] for row in cdr3map_rows}

  allrows = readCsv(enclone_file)

  cts = readCsv(celltype_file)
  for ct in ["CD4", "CD8"]:
    cells = set([row[''] for row in cts if ct in row['cell_type']])
    logger.info(f"There are {len(cells)} {ct} cells in gene expression data")

    rows = [row for row in allrows if row['barcode'] in cells]
    logger.info(f"There are {len(rows)} {ct} cells in clonotype data")

    group_clono_map = {}
    total_cells = len(rows)
    for row in rows:
      group_id = int(row['group_id'])
      clono_map = group_clono_map.setdefault(group_id, {})
      clonotype_id = int(row['clonotype_id'])
      clono_map.setdefault(clonotype_id, []).append(row)

    clono_types = []
    samples = []
    for group_id in sorted(group_clono_map.keys()):
      clono_map = group_clono_map[group_id]
      for clonotype_id in sorted(clono_map.keys()):
        cid = f"g{str(group_id).zfill(5)}c{str(clonotype_id).zfill(2)}"
        rows = clono_map[clonotype_id]
        frequency = len(rows)
        proportion = frequency / total_cells

        add_clono_data(cid, frequency, proportion, rows, cdr3map, samples, clono_types)

    samples = sorted([s for s in set(samples)])
    print(samples)

    assert total_cells == sum(clono_data['frequency'] for clono_data in clono_types)

    write_to_file(output_file_prefix + "." + ct + ".clonotype.csv", clono_types, samples)
    
  logger.info("done")

def main():
  parser = argparse.ArgumentParser(description="enclone to clonotypes",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = True
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help="Input enclone file", required=NOT_DEBUG)
  parser.add_argument('-c', '--chain', action='store', nargs='?', help="Input cdr3 chain file", required=NOT_DEBUG)
  parser.add_argument('--celltype_file', action='store', nargs='?', help="Input cell type file", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output clonotype file", required=NOT_DEBUG)
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  
  if DEBUG:
    args.input="/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_02_enclone/result/AG3669.pchain4.pcell.csv"
    args.chain="/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_01_merge/result/all_contig_annotations.json.cdr3"
    args.celltype_file="/data/h_gelbard_lab/projects/20220508_scRNA_3669/seurat_merge_03_choose_res/result/AG3669.meta.csv"
    args.output="/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_03_convert/result/AG3669"

  check_file(args.input, parser)
  check_file(args.chain, parser)
  check_file(args.celltype_file, parser)

  logger = initialize_logger()
  convert(logger, args.input, args.chain, args.celltype_file, args.output)
  
if __name__ == "__main__":
    main()
