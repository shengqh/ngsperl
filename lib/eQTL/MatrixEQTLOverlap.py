import argparse
import sys
import logging
import os
import re
import csv
from PlinkUtils import PlinkSNPItem, readPlinkSNP
from MatrixEQTLUtils import MatrixEQTLItem, readMatrixEQTLResult, fillMatrixEQTL
from GTExUtils import GTExItem, readGTExResult
from os import listdir
from os.path import isfile, join

def checkFileExist(fileName, logger):
  if not os.path.isfile(fileName):
    logger.error("file not exists: " + fileName)
    sys.exit()

def findOverlap(inputFile, inputBimFile, inputName, compareFile, compareBimFile, compareName, outputFile, logger):
  checkFileExist(inputFile, logger)
  checkFileExist(inputBimFile, logger)
  checkFileExist(compareFile, logger)
  checkFileExist(compareBimFile, logger)
  
  inputResult = readMatrixEQTLResult(inputFile)
  fillMatrixEQTL(inputResult, inputBimFile)
  logger.info("input eQTL = %d" % len(inputResult))

  compareResult = readMatrixEQTLResult(compareFile)
  fillMatrixEQTL(compareResult, compareBimFile)
  logger.info("compare eQTL = %d" % len(compareResult))

  pvalues = [0.01, 0.001, 0.0001, 0.00001, 0.000001, -0.05]
  tmpFile = outputFile + ".tmp"
  with open(tmpFile, "w") as f:
    f.write("SNP\tLocus\tGene\t%s_Major\t%s_Minor\t%s_Pvalue\t%s_FDR\t%s_Beta\t%s_Major\t%s_Minor\t%s_Pvalue\t%s_FDR\t%s_Beta\tBetaMatch\n" % (inputName, inputName, inputName, inputName, inputName, compareName, compareName, compareName, compareName, compareName));
    
    with open(outputFile + ".info", "w") as fi:
      fi.write("pValue\t%s\t%s\tOverlap\tOverlapRate\tMatched\tConflicting_eQTL\teQTL_inconsistency\n" % (inputName, compareName))
      
      for pvalue in pvalues:
        matched = 0
        unmatched = 0
        if pvalue < 0:
          fdr = -pvalue 
          logger.info("fdr = %f" % fdr)
          curInputResultMap = {p.SNP + ":" + p.Gene:p for p in inputResult if p.FDR <= fdr}
          curCompareResultMap = {p.SNP + ":" + p.Gene:p for p in compareResult if p.FDR <= fdr}
        else:
          logger.info("pvalue = %f" % pvalue)
          curInputResultMap = {p.SNP + ":" + p.Gene:p for p in inputResult if p.Pvalue <= pvalue}
          curCompareResultMap = {p.SNP + ":" + p.Gene:p for p in compareResult if p.Pvalue <= pvalue}
          
        keys = sorted(set(curInputResultMap.keys()).union(curCompareResultMap.keys()))
      
        for iKey in keys:
          if iKey in curCompareResultMap and iKey in curInputResultMap:
            i = curInputResultMap[iKey]
            c = curCompareResultMap[iKey]
            betaMatch = False;
            if (i.MajorAllele == c.MajorAllele) and (i.MinorAllele == c.MinorAllele):
              betaMatch = i.Beta * c.Beta > 0
            elif (i.MajorAllele == c.MinorAllele) and (i.MinorAllele == c.MajorAllele):
              betaMatch = i.Beta * c.Beta < 0
            else:
              betaMatch = "Unknown";
            
            if betaMatch:
              matched = matched + 1
            else:
              unmatched = unmatched + 1
            
            if pvalue == pvalues[0]:
              f.write("%s\t%s\t%s\t%s\t%s\t%.2E\t%.2E\t%f\t%s\t%s\t%.2E\t%.2E\t%f\t%r\n" % (
                i.SNP,
                i.Locus,
                i.Gene,
                i.MajorAllele,
                i.MinorAllele,
                i.Pvalue,
                i.FDR,
                i.Beta,
                c.MajorAllele,
                c.MinorAllele,
                c.Pvalue,
                c.FDR,
                c.Beta,
                betaMatch));
          elif iKey in curInputResultMap:
            i = curInputResultMap[iKey]
            if pvalue == pvalues[0]:
              f.write("%s\t%s\t%s\t%s\t%s\t%.2E\t%.2E\t%f\t\t\t\t\t\t\n" % (
                i.SNP,
                i.Locus,
                i.Gene,
                i.MajorAllele,
                i.MinorAllele,
                i.Pvalue,
                i.FDR,
                i.Beta));
          else:
            c = curCompareResultMap[iKey]
            if pvalue == pvalues[0]:
              f.write("%s\t%s\t%s\t\t\t\t\t\t%s\t%s\t%.2E\t%.2E\t%f\t\n" % (
                c.SNP,
                c.Locus,
                c.Gene,
                c.MajorAllele,
                c.MinorAllele,
                c.Pvalue,
                c.FDR,
                c.Beta));
          
        overlap = matched + unmatched
        minEqtl = min(len(curInputResultMap), len(curCompareResultMap))
        
        if minEqtl == 0:
          fi.write("%f\t%d\t%d\t%d\t%f\t%d\t%d\t%f\n" % (pvalue, len(curInputResultMap), len(curCompareResultMap),
                                                                   0,
                                                                   0,
                                                                   0,
                                                                   0,
                                                                   0))
        else:
          overlapRate = 1.0 * overlap / min(len(curInputResultMap), len(curCompareResultMap))
          fi.write("%f\t%d\t%d\t%d\t%f\t%d\t%d\t%f\n" % (pvalue, len(curInputResultMap), len(curCompareResultMap),
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
  parser = argparse.ArgumentParser(description="find overlap between MatrixEQTL results",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='MatrixEQTL file', required=NOT_DEBUG)
  parser.add_argument('-b', '--input_bim_file', action='store', nargs='?', help='Plink bim file', required=NOT_DEBUG)
  parser.add_argument('--input_name', action='store', nargs='?', help='Input name', required=NOT_DEBUG)
  parser.add_argument('-c', '--compare_file', action='store', nargs='?', help='Compare MatrixEQTL file', required=NOT_DEBUG)
  parser.add_argument('-p', '--compare_bim_file', action='store', nargs='?', help='Compare plink bim file', required=NOT_DEBUG)
  parser.add_argument('--compare_name', action='store', nargs='?', help='Compare name', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help='Output file', required=NOT_DEBUG)

  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/cqs/shengq2/guoyan/matrixEqtlNormalPaired/result/normal_paired_BRCA.cis.txt"
    args.input_bim_file = "/scratch/cqs/shengq2/guoyan/dataPreparationNormalPaired/result/BRCA/normal_paired_BRCA_snp.bim"
    args.input_name = "NORMAL"
    args.compare_file = "/scratch/cqs/shengq2/guoyan/matrixEqtlTumorPaired/result/tumor_paired_BRCA.cis.txt"
    args.compare_bim_file = "/scratch/cqs/shengq2/guoyan/dataPreparationTumorPaired/result/BRCA/tumor_paired_BRCA_snp.bim"
    args.compare_name = "TUMOR"
    args.output = "/scratch/cqs/shengq2/guoyan/normal_tumor_BRCA.comparison.tsv"
  
  logger = logging.getLogger('MatrixEQTL_Overlap')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  findOverlap(args.input, args.input_bim_file, args.input_name, args.compare_file, args.compare_bim_file, args.compare_name, args.output, logger)
  
if __name__ == "__main__":
    main()
