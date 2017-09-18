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

def findOverlap(inputFile, inputBimFile, inputName, compareFile, compareBimFile, compareName, outputFile, logger):
  inputMap = {p.Name:p for p in readPlinkSNP(inputBimFile) if p.Chromosome != "0"}
  inputResult = readMatrixEQTLResult(inputFile)
  logger.info("input eQTL = %d" % len(inputResult))

  compareMap = {p.Name:p for p in readPlinkSNP(compareBimFile) if p.Chromosome != "0"}
  compareResultMap = {p.SNP + ":" + p.Gene:p for p in readMatrixEQTLResult(compareFile)}
  logger.info("compare eQTL = %d" % len(compareResultMap))

  matched = 0
  unmatched = 0
  tmpFile = outputFile + ".tmp"
  with open(tmpFile, "w") as f:
    f.write("Locus\tGene\t%s_Major\t%s_Minor\t%s_Beta\t%s_Major\t%s_Minor\t%s_Beta\tBetaMatch\n" %(inputName, inputName, inputName, compareName, compareName, compareName));
    for i in inputResult:
      iAllele = inputMap[i.SNP];
      iKey = i.SNP + ":" + i.Gene;
      if iKey in compareResultMap:
        c = compareResultMap[iKey]
        cAllele = compareMap[i.SNP]
        betaMatch = False;
        if (iAllele.MajorAllele == cAllele.MajorAllele) and (iAllele.MinorAllele == cAllele.MinorAllele):
          betaMatch = i.Beta * c.Beta > 0
        elif (iAllele.MajorAllele == cAllele.MinorAllele) and (iAllele.MinorAllele == cAllele.MajorAllele):
          betaMatch = i.Beta * c.Beta < 0
        else:
          continue;
        
        if betaMatch:
          matched = matched + 1
        else:
          unmatched = unmatched + 1
        
        f.write("%s\t%s\t%s\t%s\t%f\t%s\t%s\t%f\t%r\n" % (
            iAllele.Locus,
            i.Gene,
            iAllele.MajorAllele,
            iAllele.MinorAllele,
            i.Beta,
            cAllele.MajorAllele,
            cAllele.MinorAllele,
            c.Beta,
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
  
  findOverlap(args.input, args.input_bim_file, args.input_name, args.compare_file, args.compare_bim_file,args.compare_name, args.output, logger)
  
if __name__ == "__main__":
    main()
