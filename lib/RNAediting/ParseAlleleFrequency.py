import pysam
import argparse
import sys
import logging
import os

from itertools import groupby
from Bio import SeqIO
from asyncore import read

DEBUG=False
NOT_DEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Parsing allele frequency from BAM files.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help="Input BAM files, separated by ','")
parser.add_argument('-o', '--output', action='store', nargs='?', required=NOT_DEBUG, help="Output file")
parser.add_argument('-f', '--fasta', action='store', nargs='?', required=NOT_DEBUG, help="Sequence fasta file, should contain only one sequence")
parser.add_argument('-b', '--min_base_quality', action='store', nargs='?', type=int, default=20, help="Minimum base quality")

args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/cqs/shengq2/emeson/20170831_rnaediting/bowtie/result/HT3.bam"
  args.fasta = "/scratch/cqs/shengq2/emeson/20170831_rnaediting/fasta/mGlu4.fasta"
  args.output = "/scratch/cqs/shengq2/emeson/20170831_rnaediting/allele/HT3.allele.tsv"

logger = logging.getLogger('parseAlleleFrequency')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

print("input=%s" % args.input)
print("fasta=%s" % args.fasta)
print("output=%s" % args.output)
print("min_base_quality=%d" % args.min_base_quality)

fastaSequences = [fs for fs in SeqIO.parse(open(args.fasta),'fasta')]
if len(fastaSequences) == 0:
  logger.error("No sequence read from file %s, exit" % args.fasta)
  sys.exit(0)
  
fastaSequence = fastaSequences[0]

basesCounts = {}

validBases = set(['A', 'T', 'G', 'C'])

readsFile = os.path.splitext(args.output)[0] + ".reads.tsv"

inputfiles = args.input.split(',')
filenames=[]
with open(readsFile, 'w') as output:
  output.write("File\tTotal\tReadFailed\tReverseStrandFailed\tValid\n")
  for inputfile in inputfiles:
    if inputfile.endswith(".bam"):
      openmode = "rb"
    else:
      openmode = "r"
      
    filename = os.path.splitext(os.path.basename(inputfile))[0]
    filenames.append(filename)
    
    basesFileCounts = {}
    for (idx, base) in enumerate(fastaSequence.seq):
      basesFileCounts[idx] = {'R':base, 'A':0, 'T':0, 'G':0, 'C':0, 'N':0}

    basesCounts[filename] = basesFileCounts
  
    logger.info("processing %s ..." % filename)
    processed = 0
    readFailedCount=0
    reversedCount=0
    samfile = pysam.Samfile(inputfile, openmode)
    try:
      for sread in samfile.fetch(until_eof=True):
        processed += 1
        if processed % 100000 == 0:
          logger.info("processed %d" % processed)
          
#        if DEBUG:
#          if processed == 1000:
#            break
#   
        if sread.is_unmapped or sread.is_secondary or sread.is_qcfail or sread.is_duplicate or sread.is_supplementary:
          readFailedCount += 1
          continue;
        
        start = sread.reference_start
        if sread.is_reverse:
          reversedCount += 1
          continue
        
        for idx, value in basesFileCounts.items():
          position = idx-start
          if position < 0:
            continue
          
          if position >= len(sread.seq):
            break
          
          if sread.query_qualities[position] < args.min_base_quality:
            continue
          
          curbase = sread.seq[position]
          if curbase in validBases:
            value['N'] = value['N']+1
            value[curbase] = value[curbase]+1
        
      logger.info("total processed %d, read failed %d, reversed %d" % (processed, readFailedCount, reversedCount ))
      output.write("%s\t%d\t%d\t%d\t%d\n" %(filename, processed, readFailedCount, reversedCount, processed - readFailedCount - reversedCount))
    finally:
      samfile.close()
      
tmpfile = args.output + ".tmp"
with open(tmpfile, 'w') as output:
  output.write("File\tPosition\tRefAllele\tRefAllelePercentage\tAltAllele\tAltAllelePercentage\tNumberOfAllele\tNumberOfA\tPercentageOfA\tNumberOfT\tPercentageOfT\tNumberOfG\tPercentageOfG\tNumberOfC\tPercentageOfC\n")
  for filename in filenames:
    basesFileCounts = basesCounts[filename]
    sortedCounts = sorted(basesFileCounts.items())
    for idx, counts in sortedCounts:
      if counts['N'] == 0:
        continue
      minor_alleles = [(allele, counts[allele] * 100.0 / counts['N']) for allele in ['A', 'T', 'G', 'C'] if allele != counts['R']]
      minor_alleles_sorted = sorted(minor_alleles, key=lambda allele:allele[1])
      minor_allele = minor_alleles_sorted[len(minor_alleles_sorted) - 1]
      for allele in ['A', 'T', 'G', 'C']:
        if allele != counts['R']:
          output.write("%s\t%d\t%s\t%.2f\t%s\t%.2f\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\n" %(filename, idx + 1, counts['R'], counts[counts['R']] * 100.0 / counts['N'], minor_allele[0], minor_allele[1], counts['N'], counts['A'], counts['A'] * 100.0 / counts['N'], counts['T'], counts['T'] * 100.0 / counts['N'], counts['G'], counts['G'] * 100.0 / counts['N'],counts['C'], counts['C'] * 100.0 / counts['N']))
          break
  if os.path.isfile(args.output):
    os.remove(args.output)
  os.rename(tmpfile, args.output)

logger.info("done.")
