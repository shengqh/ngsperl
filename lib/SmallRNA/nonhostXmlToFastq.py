import sys
import os
import logging
import argparse
import gzip
from Bio import SeqIO
from CountXmlUtils import readCountXmlQueryNames

DEBUG = 0

if DEBUG:
  inputFile="/scratch/cqs/shengq2/vickers/20180410_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_genome/bowtie1_bacteria_group1_pm_count/result/Urine_WT_14/Urine_WT_14.bam.count.mapped.xml"
  fastqFile = "/scratch/cqs/shengq2/vickers/20180410_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_unmapped_reads/result/Urine_WT_14_clipped_identical.unmapped.fastq.gz"
  outputFile="/scratch/cqs/shengq2/temp/bacteria_group1_Urine_WT_14.fastq.gz"
else:
  parser = argparse.ArgumentParser(description="Extract mapped reads to Fastq.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count xml file')
  parser.add_argument('-f', '--fastq', action='store', nargs='?', help="Original fastq file")
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output fastq file")

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  fastqFile = args.fastq
  outputFile = args.output

logger = logging.getLogger('xmlToFastq')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("Reading queries from " + inputFile)
mapped = readCountXmlQueryNames(inputFile)

with gzip.open(fastqFile, "rt") as inHandle:
  fastq_parser = SeqIO.parse(inHandle, "fastq") 
  wanted = (rec for rec in fastq_parser if rec.id in mapped)
  
  logger.info("Writing fastq to " + outputFile)
  with gzip.open(outputFile, "wt") as outHandle:
    SeqIO.write(wanted, outHandle, "fastq")

