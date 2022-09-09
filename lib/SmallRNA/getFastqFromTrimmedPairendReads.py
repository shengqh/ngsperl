import sys
import gzip
import os
import logging
import argparse
import seaborn as sns
from Bio.Seq import Seq
from difflib import SequenceMatcher

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def pass_length_restriction(minReadLength, maxReadLength, seq):
  lseq = len(seq)
  return lseq >= minReadLength and lseq <= maxReadLength

def check(logger, read1, read2, minReadLength, maxReadLength, minSimilarityRatio, output):
  sns.set_theme(style="whitegrid")

  ratios = []
  with open(output + ".txt", "wt") as fs1:
    fs1.write("read1\tread2\trc_read2\tsimilarityRatio\n")
    with gzip.open(output, "wt") as fout:
      read_count = 0
      with gzip.open(read1, "rt") as fin1:
        with gzip.open(read2, "rt") as fin2:
          while(True):
            qname1 = fin1.readline()
            qname2 = fin2.readline()

            if not qname1:
              break

            read_count += 1
            if read_count % 1000000 == 0:
              logger.info("%s" % read_count)

            seq1 = fin1.readline().rstrip()
            seq2 = fin2.readline().rstrip()

            m1 = fin1.readline()
            m2 = fin2.readline()

            s1 = fin1.readline()
            s2 = fin2.readline()

            if pass_length_restriction(minReadLength, maxReadLength, seq1) and pass_length_restriction(minReadLength, maxReadLength, seq2):
              seq = Seq(seq2)
              rc_seq = str(seq.reverse_complement())
              ss = similar(seq1, rc_seq)
              fs1.write(f"{seq1}\t{seq2}\t{rc_seq}\t{ss}\n")
              ratios.append(ss)
              if ss >= minSimilarityRatio:
                if len(rc_seq) <= len(seq1):
                  fout.write(f"{qname1}")
                  fout.write(f"{seq1}\n")
                  fout.write(f"{m1}")
                  fout.write(f"{s1}")
                else:
                  fout.write(f"{qname2}")
                  fout.write(f"{rc_seq}\n")
                  fout.write(f"{m2}")
                  fout.write(f"{s2}")
  ax = sns.violinplot(y=ratios)
  ax.figure.savefig(output + ".png", dpi=300)

DEBUG = False

if DEBUG:
  read1="/scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X6129CP9_clipped.1.fastq.short.gz"
  read2="/scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X6129CP9_clipped.2.fastq.short.gz"
  minReadLength=16
  maxReadLength=50
  minSimilarityRatio=0.7
  output="/scratch/weissvl/shengq2/20220819_rnaseq_1508_smallRNA/fastq/X6129CP9.mirna.fastq.gz"
else:
  parser = argparse.ArgumentParser(description="Generate smallRNA NTA read for Fastq file.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('--read1', action='store', nargs='?', help='Input Fastq read1 file')
  parser.add_argument('--read2', action='store', nargs='?', help='Input Fastq read2 file')
  parser.add_argument('--minReadLength', action='store', nargs='?', type=int, default=16, help="Minimum read length")
  parser.add_argument('--maxReadLength', action='store', nargs='?', type=int, default=50, help="Maximum read length")
  parser.add_argument('--minSimilarityRatio', action='store', nargs='?', type=float, default=0.7, help="Minimum similarity ratio between read1 and read2")
  parser.add_argument('-o', '--output', action='store', nargs='?', help='Output file')

  args = parser.parse_args()
  
  print(args)
  
  read1 = args.read1
  read2 = args.read2
  minReadLength = args.minReadLength
  maxReadLength = args.maxReadLength
  minSimilarityRatio = args.minSimilarityRatio
  output = args.output

logger = logging.getLogger('extract')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

check(logger, read1, read2, minReadLength, maxReadLength, minSimilarityRatio, output)