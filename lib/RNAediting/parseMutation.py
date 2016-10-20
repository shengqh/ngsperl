import pysam
import argparse
import sys
import logging
import os

from itertools import groupby
from Bio import SeqIO
from asyncore import read

parser = argparse.ArgumentParser(description="Parsing mismatch combination in BAM files.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', required=True, help="Input BAM files, separated by ','")
parser.add_argument('-o', '--output', action='store', nargs='?', required=True, help="Output file")
parser.add_argument('-f', '--fasta', action='store', nargs='?', required=True, help="Sequence fasta file, should contain only one sequence")
parser.add_argument('-p', '--positions', action='store', nargs='?', required=True, help="Position of interesting, 1-based and separated by ','")

args = parser.parse_args()

# args.input = "/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Cerebellum-Rat1_S01.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Cerebellum-Rat2_S11.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Cerebellum-Rat3_S21.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Colon-Rat1_S07.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Colon-Rat2_S17.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Colon-Rat3_S27.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Cortex-Rat1_S02.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Cortex-Rat2_S12.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Cortex-Rat3_S22.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Hippocampus-Rat1_S03.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Hippocampus-Rat2_S13.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Hippocampus-Rat3_S23.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Hypothalamus-Rat1_S04.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Hypothalamus-Rat2_S14.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Hypothalamus-Rat3_S24.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Kidney-Rat1_S09.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Kidney-Rat2_S19.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Kidney-Rat3_S29.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Lung-Rat1_S06.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Lung-Rat2_S16.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Lung-Rat3_S26.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Pancreas-Rat1_S10.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Pancreas-Rat2_S20.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Pancreas-Rat3_S30.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Stomach-Rat1_S08.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Stomach-Rat2_S18.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Stomach-Rat3_S28.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Striatum-Rat1_S05.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Striatum-Rat2_S15.bam,/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/Striatum-Rat3_S25.bam";
# args.fasta = "/workspace/shengq1/guoyan/20160922_rnaediting/database/rat_GRM4.fasta"
# args.positions = "72, 75, 231, 246, 331"

positions=[int(ps.strip()) for ps in args.positions.split(',')]

logger = logging.getLogger('parseMutation')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

fastaSequences = [fs for fs in SeqIO.parse(open(args.fasta),'fasta')]
refBases = [fastaSequences[0].seq[idx-1] for idx in positions]
refString = "".join(refBases)

#args.output = "/workspace/shengq1/guoyan/20160922_rnaediting/rnaediting.mutation." + refString + ".tsv"
summaryFile = os.path.splitext(args.output)[0] + ".summary.tsv"
readsFile = os.path.splitext(args.output)[0] + ".reads.tsv"

mutations = []
inputfiles = args.input.split(',')

with open(readsFile, 'w') as output:
  output.write("File\tTotal\tFailed\tReverseStrand\tProcessed\n")
  for inputfile in inputfiles:
    if inputfile.endswith(".bam"):
      openmode = "rb"
    else:
      openmode = "r"
      
    filename = "STDIN" if inputfile == "-" else os.path.splitext(os.path.basename(inputfile))[0]
  
    logger.info("processing %s ..." % filename)
    processed = 0
    failedCount=0
    reversedCount=0
    samfile = pysam.Samfile(inputfile, openmode)
    try:
      mutationMap = {}
      for read in samfile.fetch(until_eof=True):
        processed += 1
        if processed % 1000000 == 0:
          logger.info("processed %d" % processed)
  
        if read.is_unmapped or read.is_secondary or read.is_qcfail or read.is_duplicate or read.is_supplementary:
          failedCount += 1
          continue;
        
        start = read.reference_start
        if read.is_reverse:
          reversedCount += 1
          continue
    
        bases = [read.seq[idx-start-1] for idx in positions]
        baseString = "".join(bases)
        
        if(mutationMap.has_key(baseString)):
          mutationMap[baseString] = mutationMap[baseString] + 1
        else:
          mutationMap[baseString] = 1
        
      logger.info("total processed %d, failed %d, reversed %d" % (processed, failedCount, reversedCount))
      output.write("%s\t%d\t%d\t%d\t%d\n" %(filename, processed, failedCount, reversedCount, processed - failedCount - reversedCount))
      
      countSum = sum(mutationMap.values())
      for key in sorted(mutationMap, key=mutationMap.get, reverse=True):
        perc = mutationMap[key] * 100.0 / countSum
        if perc > 0.01:
          mutations.append([filename, key, mutationMap[key], perc])
    finally:
      samfile.close()
      
if args.output == "-":
  output = sys.stdout
else:
  tmpfile = args.output + ".tmp"
  output = open(tmpfile, 'w')

try:
  output.write("File\tMutation\tReadCount\tPercentage\n")
  for mutation in mutations:
    output.write("%s\t%s\t%d\t%.2f\n" % (mutation[0], mutation[1], mutation[2], mutation[3]))
    
  if args.output != "-":
    if os.path.isfile(args.output):
      os.remove(args.output)
    os.rename(tmpfile, args.output)
finally:
  output.close()

groups = []
mutations.sort(key=lambda mutation: mutation[1])
for k, g in groupby(mutations, lambda mutation: mutation[1]):
  groups.append([k, list(g)])
  
samples = sorted(set(mutation[0] for mutation in mutations))

groups.sort(key=lambda g:sum(x[3] for x in g[1]), reverse=True)

with open(summaryFile, 'w') as output:
  output.write("Sequence\tIsReference\t%s\n" % ("\t".join(samples)))
  for g in groups:
    output.write(g[0])
    if g[0] == refString :
      output.write("\tTrue")
    else:
      output.write("\tFalse")
    for sample in samples:
      mus = [mu for mu in g[1] if mu[0] == sample]
      if len(mus) == 0:
        output.write("\t0.0")
      else:
        output.write("\t%.2f" % mus[0][3])
    output.write("\n")
