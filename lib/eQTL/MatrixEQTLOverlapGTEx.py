import argparse
import sys
import logging
import os
import re
import csv
from MatrixEQTLUtils import MatrixEQTLItem, readMatrixEQTLResult, fillMatrixEQTL
from GTExUtils import GTExItem, readGTExResult
from os import listdir
from os.path import isfile, join

def checkFileExist(fileName, logger):
  if not os.path.isfile(fileName):
    logger.error("file not exists: " + fileName)
    sys.exit()

def findOverlap(inputFile, bimFile, gtexFile, outputFile, logger):
  checkFileExist(inputFile, logger)
  checkFileExist(bimFile, logger)
  checkFileExist(gtexFile, logger)

  gtex = readGTExResult(gtexFile)
  gtex = [g for g in gtex if len(g.RefAllele) == 1 and len(g.AltAllele) == 1]
  logger.info("GTEx = %d" % len(gtex))
  
  gtexMap = {};
  for g in gtex:
    if g.Key not in gtexMap:
      gtexMap[g.Key] = g;

  mResult = readMatrixEQTLResult(inputFile)
  fillMatrixEQTL(mResult, bimFile)
  logger.info("eQTL = %d" % len(mResult))

  pvalues = [0.01, 0.001, 0.0001, 0.00001, 0.000001, -0.05]
  tmpFile = outputFile + ".tmp"
  with open(tmpFile, "w") as f:
    f.write("Locus\tGene\tTCGA_Major\tTCGA_Minor\tTCGA_Pvalue\tTCGA_FDR\tTCGA_Beta\tGTex_Ref\tGTex_Alt\tGTex_Pvalue\tGTex_FDR\tGTex_Beta\tBetaMatch\n");
    with open(outputFile + ".info", "w") as fi:
      fi.write("pValue\tMatrixEQTL\tGTEx\tOverlap\tOverlapRate\tMatched\tConflicting_eQTL\teQTL_inconsistency\n")
      for pvalue in pvalues:
        matched = 0
        unmatched = 0
        if pvalue < 0:
          fdr = -pvalue 
          logger.info("fdr = %f" % fdr)
          curInputResult = [p for p in mResult if p.FDR <= fdr]
          curCompareResultMap = {g:v for g,v in gtexMap.iteritems() if v.FDR <= fdr}
        else:
          logger.info("pvalue = %f" % pvalue)
          matched = 0
          unmatched = 0
          curInputResult = [p for p in mResult if p.Pvalue <= pvalue]
          curCompareResultMap = {g:v for g,v in gtexMap.iteritems() if v.Pvalue <= pvalue}
      
        for b in curInputResult:
          key = b.Locus + ":" + b.Gene;
          if key in curCompareResultMap:
            gtexItem = curCompareResultMap[key]
            betaMatch = False;
            if (b.MajorAllele == gtexItem.RefAllele) and (b.MinorAllele == gtexItem.AltAllele):
              betaMatch = b.Beta * gtexItem.Beta > 0
            elif (b.MajorAllele == gtexItem.AltAllele) and (b.MinorAllele == gtexItem.RefAllele):
              betaMatch = b.Beta * gtexItem.Beta < 0
            else:
              continue;
            
            if betaMatch:
              matched = matched + 1
            else:
              unmatched = unmatched + 1
            
            if pvalue == pvalues[0]: 
              f.write("%s\t%s\t%s\t%s\t%f\t%.2E\t%s\t%s\t%f\t%.2E\t%r\n" %(
                  b.Locus,
                  b.Gene,
                  b.MajorAllele,
                  b.MinorAllele,
                  b.Beta,
                  b.Pvalue,
                  gtexItem.RefAllele,
                  gtexItem.AltAllele,
                  gtexItem.Beta,
                  gtexItem.Pvalue,
                  betaMatch));
                  
        overlap = matched + unmatched
        minEqtl = min(len(curInputResult), len(curCompareResultMap))
        
        if minEqtl == 0:
          fi.write("%f\t%d\t%d\t%d\t%f\t%d\t%d\t%f\n" % (pvalue, len(curInputResult), len(curCompareResultMap),
                                                                   0,
                                                                   0,
                                                                   0,
                                                                   0,
                                                                   0))
        else:
          overlapRate = 1.0 * overlap / min(len(curInputResult), len(curCompareResultMap))
          fi.write("%f\t%d\t%d\t%d\t%f\t%d\t%d\t%f\n" % (pvalue, len(curInputResult), len(curCompareResultMap),
                                                                   overlap,
                                                                   overlapRate,
                                                                   matched,
                                                                   unmatched,
                                                                   (unmatched * 1.0 / (matched + unmatched))))
            
  if os.path.isfile(outputFile):
    os.remove(outputFile)
  os.rename(tmpFile, outputFile)
  
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
    args.input = "/scratch/cqs/shengq2/guoyan/matrixEqtl_rnaTP_snpTP/result/rnaTP_snpTP_BRCA.cis.txt"
    args.bim_file = "/scratch/cqs/shengq2/guoyan/dataPreparation_rnaTP_snpTP/result/BRCA/rnaTP_snpTP_BRCA_snp.bim"
    args.gtex_file = "/scratch/cqs/shengq2/references/GTEx/Breast_Mammary_Tissue_Analysis.snpgenes"
    args.output = "/scratch/cqs/shengq2/guoyan/testOverlap.tsv"
  
  logger = logging.getLogger('MatrixEQTL_GTEx')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  findOverlap(args.input, args.bim_file, args.gtex_file, args.output, logger)
  
if __name__ == "__main__":
    main()
