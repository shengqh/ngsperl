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

def readCsv(filename, sep=","):
  allrows = []
  with open(filename, "rt") as fin:
    reader = csv.DictReader(fin, delimiter=sep)
    for row in reader:
      allrows.append(row)
  return(allrows)

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

def has_any_celltype(ct_set, cell_type):
  for ct in ct_set:
    if ct in cell_type:
      return(True)
  return(False)

def convert(logger, enclone_file, chain_file, celltype_file, output_file):
  cdr3map_rows = readCsv(chain_file, sep="\t")
  cdr3map = {row['cdr3']:row['chain'] for row in cdr3map_rows}

  allrows = readCsv(enclone_file)

  ct_map = {'all': set()}
  if celltype_file != None:
    cts = readCsv(celltype_file)

    tcell=set(row['cell_type'] for row in cts if has_any_celltype(['CD4', 'CD8', 'T cells'], row['cell_type']))
    if len(tcell) > 0:
      ct_map["T_cells"] = tcell

      cd4=set(row['cell_type'] for row in cts if has_any_celltype(['CD4'], row['cell_type']))
      if len(cd4) > 0:
        ct_map["CD4"] = cd4
        if len(cd4) > 1:
          for cd4sub in cd4:
            ct_map[cd4sub.replace(" ", "_")] = set([cd4sub])

      cd8=set(row['cell_type'] for row in cts if has_any_celltype(['CD8'], row['cell_type']))
      if len(cd8) > 0:
        ct_map["CD8"] = cd8
        if len(cd8) > 1:
          for cd8sub in cd8:
            ct_map[cd8sub.replace(" ", "_")] = set([cd8sub])

  with open(output_file, "wt") as flist:
    flist.write("cell_type\tmeta_file\n")
    for ct_name in ct_map.keys():
      if ct_name != 'all':
        ct_set = ct_map[ct_name]
        cells = set([row[''] for row in cts if has_any_celltype(ct_set, row['cell_type'])])
        logger.info(f"There are {len(cells)} {ct_name} cells in gene expression data")

        rows = [row for row in allrows if row['barcode'] in cells]
        logger.info(f"There are {len(rows)} {ct_name} cells in clonotype data")
      else:
        rows = allrows

      clone_sample_cell_map = {}
      for row in rows:
        group_id = int(row['group_id'])

        barcode = row['barcode']
        m = re.search('(.+)_', barcode)
        sample = m.group(1)

        clone_sample_cell_map.setdefault(sample, {}).setdefault(group_id, []).append(row)

      samples = sorted([s for s in clone_sample_cell_map.keys()])
      print(samples)

      output_folder = os.path.dirname(os.path.abspath(output_file))
      sub_folder = os.path.join(output_folder, ct_name)
      if not os.path.exists(sub_folder):
        os.mkdir(sub_folder)
      
      meta_file = os.path.join(sub_folder, ct_name + ".metadata.txt")
      flist.write(f"{ct_name}\t{meta_file}\n")

      with open(meta_file, "wt") as fmeta:
        fmeta.write("Sample\tCondition\n")
        for sample in samples:
          fname = f"{ct_name}_{sample}"
          fmeta.write(f"{fname}\t{sample}\n")
          group_map = clone_sample_cell_map[sample]
          sample_file = os.path.join(sub_folder, fname+ ".csv")
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
  parser.add_argument('--celltype_file', action='store', nargs='?', help="Input cell type file", required=False)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output meta file", required=NOT_DEBUG)
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  
  if DEBUG:
    args.input="/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_02_enclone/result/AG3669.pchain4.pcell.csv"
    args.chain="/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_01_merge/result/all_contig_annotations.json.cdr3"
    args.celltype_file="/data/h_gelbard_lab/projects/20220508_scRNA_3669/seurat_merge_03_choose_res/result/AG3669.meta.csv"
    args.output="/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_05_consensus/result/AG3669.meta.list"

  check_file(args.input, parser)
  if args.celltype_file != None:
    check_file(args.celltype_file, parser)
  check_file(args.chain, parser)

  logger = initialize_logger()
  convert(logger, args.input, args.chain, args.celltype_file, args.output)
  
if __name__ == "__main__":
    main()
