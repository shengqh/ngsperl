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

def pass_length_restriction(minReadLength, seq):
  lseq = len(seq)
  return lseq >= minReadLength

def output_read(fout, qname1, seq1, m1, s1, sequence_count): 
  fout.write(f"{qname1}")
  fout.write(f"{seq1}\n")
  fout.write(f"{m1}")
  fout.write(f"{s1}\n")
  if seq1 in sequence_count:
    sequence_count[seq1] = sequence_count[seq1] + 1
  else:
    sequence_count[seq1] = 1

def check(logger, read1, read2, minReadLength, minSimilarityRatio, output):
  sns.set_theme(style="whitegrid")

  len_not_equals = 0
  too_short = 0
  too_long = 0
  less_similarity = 0
  has_N = 0
  ambiguous = 0
  identical = 0
  r1_in_r2 = 0
  r2_in_r1 = 0
  by_identical = 0
  by_fix_N = 0
  read_count = 0

  ratios = []
  sequence_count = {}
  waiting = []
  #length_dic = {}
  with gzip.open(output + ".discarded.txt.gz", "wt") as fs2, gzip.open(output + ".txt.gz", "wt") as fs1, gzip.open(output, "wt") as fout, gzip.open(read1, "rt") as fin1, gzip.open(read2, "rt") as fin2:
    fs2.write("qname\tread1\tread2\trc_read2\tsimilarityRatio\treason\n")
    fs1.write("read1\tread2\trc_read2\tsimilarityRatio\n")
    while(True):
      qname1 = fin1.readline()
      qname2 = fin2.readline()

      if not qname1:
        break

      qname = qname1.rstrip()

      read_count += 1
      if read_count % 1000000 == 0:
        logger.info("%s" % read_count)

      seq1 = fin1.readline().rstrip()
      seq2 = fin2.readline().rstrip()

      m1 = fin1.readline()
      m2 = fin2.readline()

      s1 = fin1.readline().rstrip()
      s2 = fin2.readline().rstrip()

      seq = Seq(seq2)
      rc_seq = str(seq.reverse_complement())
      ss = similar(seq1, rc_seq)

      ratios.append(ss)
      fs1.write(f"{seq1}\t{seq2}\t{rc_seq}\t{ss}\n")

      min_len = min(len(seq1), len(seq2))
      if min_len < minReadLength:
        too_short += 1
        fs2.write(f"{qname}\t{seq1}\t{seq2}\t{rc_seq}\t{ss}\ttoo_short\n")
        continue

      if ss == 1:
        identical += 1
        output_read(fout, qname1, seq1, m1, s1, sequence_count)
        continue

      if len(seq1) != len(rc_seq):
        if (len(seq1) < len(rc_seq)) and rc_seq.endswith(seq1): #read1 correct, read2 adapter trimming failed
          r1_in_r2 += 1
          output_read(fout, qname1, seq1, m1, s1, sequence_count)
          continue

        if (len(rc_seq) < len(seq1)) and seq1.startswith(rc_seq):#read2 correct, read1 adapter trimming failed
          r2_in_r1 += 1
          output_read(fout, qname1, rc_seq, m1, s2[::-1], sequence_count)
          continue

        len_not_equals += 1
        fs2.write(f"{qname}\t{seq1}\t{seq2}\t{rc_seq}\t{ss}\tnot_equal_length\n")
        continue

      #eqaul length but not match
      if ss < minSimilarityRatio:
        less_similarity += 1
        fs2.write(f"{qname}\t{seq1}\t{seq2}\t{rc_seq}\t{ss}\tless_similarity\n")
        continue

      waiting.append([qname, seq1, m1, s1, rc_seq, s2, ss])
        
    #check waiting entris
    for v in waiting:
      qname = v[0]
      r1 = v[1]
      r2 = v[4]

      #pick the one with highest reads 
      c1 = sequence_count[r1] if r1 in sequence_count else 0
      c2 = sequence_count[r2] if r2 in sequence_count else 0

      if c1 > c2 * 2:
        by_identical += 1
        fout.write(f"{qname}\n")
        fout.write(f"{r1}\n")
        fout.write(f"{v[2]}")
        fout.write(f"{v[3]}\n")
        continue
      
      if c2 > c1 * 2:
        by_identical += 1
        fout.write(f"{qname}\n")
        fout.write(f"{r2}\n")
        fout.write(f"{v[2]}")
        fout.write(f"{v[5][::-1]}\n")
        continue

      #try remove N difference
      chars = []
      for idx in range(0, len(r1)):
        if r1[idx] == r2[idx]:
          chars.append(r1[idx])
        elif r1[idx] == 'N':
          chars.append(r2[idx])
        elif r2[idx] == 'N':
          chars.append(r1[idx])
        else:
          chars.clear()
          break

      #all N have been replaced and no other mismatch
      if len(chars) > 0:
        by_fix_N += 1
        newseq = "".join(chars)
        fout.write(f"{qname}\n")
        fout.write(f"{newseq}\n")
        fout.write(f"{v[2]}")
        fout.write(f"{v[3]}\n")
        continue

      #if there is N in either reads, discard
      if 'N' in r1 or 'N' in r2:
        has_N += 1
        fs2.write(f"{qname}\t{r1}\t{v[2]}\t{v[3]}\t{v[4]}\thas_N\n")
        continue

      ambiguous += 1
      fs2.write(f"{qname}\t{r1}\t{v[2]}\t{v[3]}\t{v[4]}\tambiguous\n")
      continue

  ax = sns.violinplot(y=ratios)
  ax.figure.savefig(output + ".png", dpi=300)

  with open(output + ".info", "wt") as fout:
    fout.write("category\treason\tcount\n")
    fout.write(f"total\ttotal\t{read_count}\n")
    fout.write(f"discard\ttoo_short\t{too_short}\n")
    fout.write(f"valid\tidentical\t{identical}\n")
    fout.write(f"valid\tr1_in_r2\t{r1_in_r2}\n")
    fout.write(f"valid\tr2_in_r1\t{r2_in_r1}\n")
    fout.write(f"discard\tnot_equal_length\t{len_not_equals}\n")
    fout.write(f"discard\tless_similarity\t{less_similarity}\n")
    fout.write(f"valid\tby_identical\t{by_identical}\n")
    fout.write(f"valid\tby_fix_N\t{by_fix_N}\n")
    fout.write(f"discard\tboth_has_N\t{has_N}\n")
    fout.write(f"discard\tambiguous\t{ambiguous}\n")

DEBUG = False

if DEBUG:
  read1="/scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X6129CP9_clipped.1.fastq.short.gz"
  read2="/scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/cutadapt/result/X6129CP9_clipped.2.fastq.short.gz"
  minReadLength=16
  minSimilarityRatio=0.7
  output="/scratch/weissvl/shengq2/20220819_rnaseq_1508_smallRNA/fastq/X6129CP9.mirna.fastq.gz"
else:
  parser = argparse.ArgumentParser(description="Generate smallRNA NTA read for Fastq file.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('--read1', action='store', nargs='?', help='Input Fastq read1 file')
  parser.add_argument('--read2', action='store', nargs='?', help='Input Fastq read2 file')
  parser.add_argument('--minReadLength', action='store', nargs='?', type=int, default=16, help="Minimum read length")
  parser.add_argument('--minSimilarityRatio', action='store', nargs='?', type=float, default=0.7, help="Minimum similarity ratio between read1 and read2")
  parser.add_argument('-o', '--output', action='store', nargs='?', help='Output file')

  args = parser.parse_args()
  
  print(args)
  
  read1 = args.read1
  read2 = args.read2
  minReadLength = args.minReadLength
  minSimilarityRatio = args.minSimilarityRatio
  output = args.output

logger = logging.getLogger('extract')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

check(logger, read1, read2, minReadLength, minSimilarityRatio, output)