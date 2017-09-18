import argparse
import sys
import logging
import os
import re
import csv
from PlinkUtils import PlinkSNPItem, readPlinkSNP
from MatrixEQTLUtils import MatrixEQTLItem, readMatrixEQTLResult
from GTExUtils import GTExItem, readGTExResult
from os import listdir
from os.path import isfile, join

def findOverlap(inputFile, bimFile, gtexFile, outputFile, logger):
  pSNPs = readPlinkSNP(bimFile)
  pSNPs = [p for p in pSNPs if p.Chromosome != "0"]
  logger.info("SNP = %d" % len(pSNPs))

  snpMap = {};
  for t in pSNPs:
    snpMap[t.Name] = t

  gtex = readGTExResult(gtexFile)
  gtex = [g for g in gtex if len(g.RefAllele) == 1 and len(g.AltAllele) == 1]
  logger.info("GTEx = %d" % len(gtex))
  
  gtexMap = {};
  for g in gtex:
    if g.Key not in gtexMap:
      gtexMap[g.Key] = g;

  mResult = readMatrixEQTLResult(inputFile)
  logger.info("eQTL = %d" % len(mResult))

  matched = 0
  unmatched = 0
  tmpFile = outputFile + ".tmp"
  with open(tmpFile, "w") as f:
    f.write("Locus\tGene\tTCGA_Major\tTCGA_Minor\tTCGA_Beta\tGTex_Ref\tGTex_Alt\tGTex_Beta\tBetaMatch\n");
    for b in mResult:
      allele = snpMap[b.SNP];
      key = allele.Locus + ":" + b.Gene;
      if key in gtexMap:
        gtexItem = gtexMap[key]
        betaMatch = False;
        if (allele.MajorAllele == gtexItem.RefAllele) and (allele.MinorAllele == gtexItem.AltAllele):
          betaMatch = b.Beta * gtexItem.Beta > 0
        elif (allele.MajorAllele == gtexItem.AltAllele) and (allele.MinorAllele == gtexItem.RefAllele):
          betaMatch = b.Beta * gtexItem.Beta < 0
        else:
          continue;
        
        if betaMatch:
          matched = matched + 1
        else:
          unmatched = unmatched + 1
        
        f.write("%s\t%s\t%s\t%s\t%f\t%s\t%s\t%f\t%r\n" %(
            allele.Locus,
            gtexItem.Gene,
            allele.MajorAllele,
            allele.MinorAllele,
            b.Beta,
            gtexItem.RefAllele,
            gtexItem.AltAllele,
            gtexItem.Beta,
            betaMatch));
            
  if os.path.isfile(outputFile):
    os.remove(outputFile)
  os.rename(tmpFile, outputFile)
  
  with open(outputFile + ".info", "w") as f:
    f.write("Type\tCount\n");
    f.write("Matched\t%d\n" % matched)
    f.write("Unmatched\t%d\n" % unmatched)
    f.write("UnmatchedPercentage\t%4.2f%%\n" % (unmatched * 100.0 / (matched + unmatched)))
  
  logger.info("Done.")
        
def main():
  parser = argparse.ArgumentParser(description="find overlap between MatrixEQTL and GTEx",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='MatrixEQTL file', required=NOT_DEBUG)
  parser.add_argument('-b', '--bim_file', action='store', nargs='?', help='Plink bim file', required=NOT_DEBUG)
  parser.add_argument('-g', '--gtex_file', action='store', nargs='?', help='GTEx file', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help='Output file', required=NOT_DEBUG)

  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/cqs/shengq2/guoyan/matrixEqtlNormal/result/normal_BRCA.cis.txt"
    args.bim_file = "/scratch/cqs/shengq2/guoyan/bim/normal_BRCA_common_filtered.bim"
    args.gtex_file = "/scratch/cqs/shengq2/references/GTEx/Breast_Mammary_Tissue_Analysis.snpgenes"
    args.output = "/scratch/cqs/shengq2/guoyan/matrixEqtlNormal/result/normal_BRCA.cis.GTEx.tsv"
  
  logger = logging.getLogger('MatrixEQTL_GTEx')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  findOverlap(args.input, args.bim_file, args.gtex_file, args.output, logger)
  
if __name__ == "__main__":
    main()
