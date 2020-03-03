import sys
import os
import os.path
import logging
import argparse
import pysam
import subprocess
from os import path

DEBUG = 1

if DEBUG:
  inputFile="/scratch/cqs/maureen_gannon_projects/20200131_rnaseq_joseph_mmul10_byTiger_old/star_featurecount/result/GLDPH.bam"
  #inputFile="/scratch/cqs/maureen_gannon_projects/20200131_rnaseq_joseph_mmul10_byTiger_old/star_featurecount/result/CTR_33215_F_Aligned.sortedByCoord.out.bam"
  fastaFile="/scratch/cqs_share/references/ensembl/Mmul_10/Macaca_mulatta.Mmul_10.dna.primary_assembly.fa"
  region="11:6685714-6689626"
  outputFile="/scratch/cqs/maureen_gannon_projects/GAPDH_consensus_sequence.txt"
else:
  parser = argparse.ArgumentParser(description="Extract consensus sequence from BAM file.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bam file')
  parser.add_argument('-f', '--fasta', action='store', nargs='?', help='Genome fasta file')
  parser.add_argument('-l', '--region', action='store', nargs='?', help="Region in genome")
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output consensus sequence file")

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  fastaFile = args.fasta
  region = args.region
  outputFile = args.output

logger = logging.getLogger('getConsensusSequence')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

bamFile = inputFile

logger.info("Mpileup BAM file " + bamFile)
tmpFile = outputFile + ".mpileup"
if not path.exists(tmpFile):
  os.system("samtools mpileup -a -r " + region + " -f " + fastaFile + " --output " + tmpFile + " " + bamFile)

gapChar = '.'
consensusSequence = ""
with open(tmpFile, "rt") as fin:
  for line in fin:
    parts = line.split('\t')
    refBase = parts[2][0]
    bases = parts[4]
    refMap = {'A':0, 'T':0, 'G':0, 'C':0, gapChar:0}
    for base in bases:
      if base == '^' or base == '$':
        continue
      
      if base == '.' or base == ',' :
        refMap[refBase] = refMap[refBase] + 1
        continue

      if base == '>' or base == '<':
        refMap[gapChar] = refMap[gapChar] + 1
        continue
      
      refMap[base] = refMap[refBase] + 1
    
    maxValue = max(refMap.values())
    if refMap[refBase] == maxValue:
      consensusSequence = consensusSequence + refBase
      continue

    if refMap[gapChar] == maxValue:
      consensusSequence = consensusSequence + refBase.lower()
      continue
    
    for base in refMap.keys():
      if refMap[base] == maxValue:
        consensusSequence = consensusSequence + base
        break

with open(outputFile, "wt") as fout:
  fout.write(consensusSequence + "\n")

logger.info("Done")
