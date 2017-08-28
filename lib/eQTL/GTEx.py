import argparse
import sys
import logging
import os
import re
import csv
from os import listdir
from os.path import isfile, join

def findGTEx(inputFile, gtexFolder, outputFile, logger):
  header = "";
  tmpfile = outputFile + ".tmp"
  
  snps = {}
  keys = []
  with open(tmpfile, 'w') as wr:
    with open(inputFile, 'r') as f:
      for line in f:
        if line.startswith("#"):
          wr.write(line)
        elif line == "":
          continue
        else:
          header = line
          break
        
      for line in f:
        parts = line.strip().split('\t')
        key = parts[0] + "_" + parts[1] + "_" + parts[3] + "_" + parts[4]
        snps[key] = {
          "line":line, "gtex":{}
        }
        keys.append(key)
    
      logger.info("total %d SNPs" % len(keys))
    
    gtexFiles = [f for f in listdir(gtexFolder) if isfile(join(gtexFolder, f))]
    gtexFiles.sort()
    tissues = []
    for gtex in gtexFiles[1:2]:
      tissue = re.sub(r"_Analysis.*", "", gtex)
      tissues.append(tissue)
      gtexfile = join(gtexFolder, gtex)
      logger.info("reading " + gtexfile)
      with open(gtexfile, 'r') as reader: 
        sr = csv.DictReader(reader, delimiter='\t')
        for row in sr:
          key1 = row["snp_chrom"] + "_" + row["snp_pos"] + "_" + row["ref"] + "_" + row["alt"];
          key2 = row["snp_chrom"] + "_" + row["snp_pos"] + "_" + row["alt"] + "_" + row["ref"];
          key = "";
          if key1 in snps:
            key = key1
          elif key2 in snps:
            key = key2
          else:
            continue
          
          if not tissue in snps[key]["gtex"]:
            snps[key]["gtex"][tissue] = [row["gene_name"]]
          else:
            snps[key]["gtex"][tissue].append(row["gene_name"])
          print("%s => %s" % (key, snps[key]["gtex"][tissue]));

    
      
      
  if os.path.isfile(outputFile):
    os.remove(outputFile)
  os.rename(tmpfile, outputFile)
        
def main():
  parser = argparse.ArgumentParser(description="find eQTL annotation in GTEx database",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = True
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='VCF file', required=NOT_DEBUG)
  parser.add_argument('-g', '--gtex_folder', action='store', nargs='?', help='GTEx folder', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help='Output file', required=NOT_DEBUG)

  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/cqs/shengq1/evan/20170823_evan_annovar/annovar/result/Highest/Highest.assoc.linear_Top500.vcf.annovar.final.tsv"
    args.gtex_folder = "/scratch/cqs/shengq1/references/GTEx"
    args.output = "/scratch/cqs/shengq1/evan/20170823_evan_annovar/annovar/result/Highest/Highest.assoc.linear_Top500.vcf.annovar.final.GTEx.tsv"
  
  logger = logging.getLogger('gtex')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  findGTEx(args.input, args.gtex_folder, args.output, logger)
  
if __name__ == "__main__":
    main()
