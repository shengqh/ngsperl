import argparse
import logging
import os
import sys
import re
import csv
from collections import OrderedDict

def initialize_logger():
  logger = logging.getLogger('enclone_to_consensus')
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

def sortByCells(val): 
    return val['n'] 

def convert(enclone_file, chain_file, output_meta_file, logger):
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

  clone_sample_cell_map = {}
  for row in rows:
    group_id = int(row['group_id'])

    barcode = row['barcode']
    m = re.search('(.+)_', barcode)
    sample = m.group(1)

    clone_sample_cell_map.setdefault(sample, {}).setdefault(group_id, []).append(row)

  samples = sorted([s for s in clone_sample_cell_map.keys()])
  print(samples)

  output_folder = os.path.dirname(output_meta_file)
  with open(output_meta_file, "wt") as fmeta:
    fmeta.write("Sample\tCondition\n")
    for sample in samples:
      fmeta.write(f"{sample}\t{sample}\n")
      group_map = clone_sample_cell_map[sample]
      sample_file = os.path.join(output_folder, sample+ ".csv")
      with open(sample_file, "wt") as fout:
        fout.write("clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis\n")

        for group_id in sorted(list(group_map.keys())):
          clonotype_id = f"clonotype{group_id}"
          cell_list = group_map[group_id]
          row = cell_list[0]
          nchain = int(row['nchains'])
          for chain_idx in range(1, nchain+1):
            cdr3 = row[f"cdr3_aa{chain_idx}"]
            if cdr3 == "":
              continue

            consensus_id = f"{clonotype_id}_consensus_{chain_idx}"
            length = 500
            v_gene = row[f"v_name{chain_idx}"] if row[f"v_name{chain_idx}"] != "" else "None"
            d_gene = row[f"d_name{chain_idx}"] if row[f"d_name{chain_idx}"] != "" else "None"
            j_gene = row[f"j_name{chain_idx}"] if row[f"j_name{chain_idx}"] != "" else "None"
            c_gene = row[f"const{chain_idx}"]
            full_length = "TRUE"
            productive = "TRUE"
            chain = cdr3map[cdr3]
            cdr3_nt = row[f"cdr3_dna{chain_idx}"]
            r_cell_key = f"r_cell{chain_idx}"
            u_cell_key = f"u_cell{chain_idx}"
            reads = sum(int(cell[r_cell_key]) if cell[r_cell_key] != "" else 0 for cell in cell_list)
            umis = sum(int(cell[u_cell_key]) if cell[u_cell_key] != "" else 0 for cell in cell_list)

            if reads > 0:
              fout.write(f"{clonotype_id},{consensus_id},{length},{chain},{v_gene},{d_gene},{j_gene},{c_gene},{full_length},{productive},{cdr3},{cdr3_nt},{reads},{umis}\n")

  logger.info("done")

def main():
  parser = argparse.ArgumentParser(description="convert enclone to consensus_annotations",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help="Input enclone file", required=NOT_DEBUG)
  parser.add_argument('-c', '--chain', action='store', nargs='?', help="Input cdr3 chain file", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output meta file", required=NOT_DEBUG)
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  
  if DEBUG:
    args.input="/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_2_enclone/result/AG3669.pcell.pchain4.csv"
    args.chain="/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_1_merge/result/all_contig_annotations.json.cdr3"
    args.output="/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_4_convert_consensus/result/metadata.txt"

  check_file(args.input, parser)
  check_file(args.chain, parser)

  logger = initialize_logger()
  convert(args.input, args.chain, args.output, logger)
  
if __name__ == "__main__":
    main()
