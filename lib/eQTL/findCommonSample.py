import argparse
import sys
import logging
import os
import re

def findSample(rnaseqFiles, rnaseqNames, famFiles, famNames, outputFilePrefix, pattern, logger):
  commonKey = set()
  bFirst = True
  rnaseqMap = {}
  for idx in range(0, len(rnaseqFiles)):
    rnaseqFile = rnaseqFiles[idx]
    rnaseqName = rnaseqNames[idx]
    
    logger.info("reading %s ..." % (rnaseqFile))
    with open(rnaseqFile, 'r') as f:
      ids = [id for id in f.readline().strip().split('\t')[2:] if re.match(pattern, id)]
      idmap = dict((re.match(pattern, v).group(1), v) for v in ids)
      rnaseqMap[rnaseqName] = idmap
      if bFirst:
        commonKey = set(idmap)
        bFirst = False
      else:
        commonKey.intersection_update(set(idmap))
      
  snpMap = {}
  for idx in range(0, len(famFiles)):
    famFile = famFiles[idx]
    famName = famNames[idx]
    
    logger.info("reading %s ..." % (famFile))
    with open(famFile, 'r') as f:
      idmap = {}
      for line in f:
        parts = line.split(' ')
        m = re.match(pattern, parts[0])
        if not m:
          raise Exception("Pattern %s not in sample %s" % (pattern, parts[0]))
        idmap[m.group(1)] = parts[0] + " " + parts[1]
      snpMap[famName] = idmap
      commonKey.intersection_update(set(idmap))
  
  logger.info("common samples %d ..." % len(commonKey))
  commonSamples = sorted(commonKey)
  for rnaseqName in rnaseqMap.keys():
    idmap = rnaseqMap[rnaseqName]
    outputFile =outputFilePrefix + ".rnaseq." + rnaseqName + ".samples" 
    logger.info("writing %s ..." % outputFile)
    with open(outputFile, 'w') as sw:
      for commonSample in commonSamples:
        sw.write(idmap[commonSample] + "\n")
     
  for famName in snpMap.keys():
    idmap = snpMap[famName]
    outputFile =outputFilePrefix + ".snp." + famName + ".samples" 
    logger.info("writing %s ..." % outputFile)
    with open(outputFile, 'w') as sw:
      for commonSample in commonSamples:
        sw.write(idmap[commonSample] + "\n")
  
  logger.info("done.")
        
def main():
  parser = argparse.ArgumentParser(description="find common samples between RNASeq and SNP data",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG

  parser.add_argument('-r', '--rnaseqFiles', action='store', nargs='?', help='RNASeq data files, delimited by ","', required=NOT_DEBUG)
  parser.add_argument('--rnaseqNames', action='store', nargs='?', help='RNASeq data names, delimited by ",", will be used in output file', required=NOT_DEBUG)
  parser.add_argument('-f', '--famFiles', action='store', nargs='?', help='Plink fam files, delimited by ","', required=NOT_DEBUG)
  parser.add_argument('--famNames', action='store', nargs='?', help='Plink fam names, delimited by ",", will be used in output file', required=NOT_DEBUG)
  parser.add_argument('-o', '--outputPrefix', action='store', nargs='?', help='Output file prefix', required=NOT_DEBUG)
  parser.add_argument('-p', '--pattern', action='store', nargs='?', help='Name pattern to match between RNAseq and plink sample name', default="(.*)", required=False)
  
  args = parser.parse_args()

  if DEBUG:
    args.rnaseqFiles = "/scratch/cqs/shengq2/guoyan/prepareRnaseq/result/BRCA/BRCA.rnaseq2.fpkm.Solid_Tissue_Normal.tsv,/scratch/cqs/shengq2/guoyan/prepareRnaseq/result/BRCA/BRCA.rnaseq2.count.Primary_Solid_Tumor.tsv"
    args.rnaseqNames = "Solid_Tissue_Normal,Primary_Solid_Tumor"
    args.famFiles = "/scratch/cqs/shengq2/guoyan/prepareSNP/result/BRCA/BRCA_Solid_Tissue_Normal.fam,/scratch/cqs/shengq2/guoyan/prepareSNP/result/BRCA/BRCA_Primary_Solid_Tumor.fam";
    args.famNames = "Solid_Tissue_Normal,Primary_Solid_Tumor"
    args.pattern = "(TCGA........)"
    args.outputPrefix = "/scratch/cqs/shengq2/guoyan/rnaNT_rnaTP_snpNT_snpTP_BRCA_common"

  print(args)
  
  logger = logging.getLogger('findCommonSample')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  findSample(args.rnaseqFiles.split(','), args.rnaseqNames.split(','), args.famFiles.split(','), args.famNames.split(','), args.outputPrefix, args.pattern, logger)
  
if __name__ == "__main__":
    main()

