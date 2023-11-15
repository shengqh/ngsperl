import sys
import gzip
import os
import logging
import argparse

DEBUG = 0

if DEBUG:
  input="/nobackup/shah_lab/shengq2/20231114_Heart_BEAT_rnaseq_hg38_umi/merge_fastq/result/POPDC2_C11.1.fastq.gz,/nobackup/shah_lab/shengq2/20231114_Heart_BEAT_rnaseq_hg38_umi/merge_fastq/result/POPDC2_C11.2.fastq.gz"
  output="/nobackup/shah_lab/shengq2/20231114_Heart_BEAT_rnaseq_hg38_umi/umitools_extract/result/POPDC2_C11.umi.1.fastq.gz,/nobackup/shah_lab/shengq2/20231114_Heart_BEAT_rnaseq_hg38_umi/umitools_extract/result/POPDC2_C11.umi.2.fastq.gz"
else:
  parser = argparse.ArgumentParser(description="Extract UMI for Pico v3 SMART UMI reads.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Fastq files, separated by comma')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output Fastq files, separated by comma")

  args = parser.parse_args()
  
  print(args)
  
  input = args.input
  output = args.output

inputFiles = input.split(",")
outputFiles = output.split(",")

read1file = inputFiles[0]
read2file = inputFiles[1]

outputFile1 = outputFiles[0]
outputFile2 = outputFiles[1]

logger = logging.getLogger('extract_umi')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info(f"Processing queries in {read1file}, {read2file} ...")

with gzip.open(read1file, 'rt') as f1, gzip.open(read2file, 'rt') as f2:
  with gzip.open(outputFile1, "wt") as fout1, gzip.open(outputFile2, 'wt') as fout2:
    readCount = 0
    while True:
      header1 = f1.readline().strip()
      header2 = f2.readline().strip()
      if '' == header1:
        break

      if not header1.startswith("@"):
        continue

      header1parts = header1.split(' ', 1)
      header2parts = header2.split(' ', 1)

      seq1 = f1.readline().strip()
      ignore1 = f1.readline().strip()
      score1 = f1.readline().strip()

      seq2 = f2.readline().strip()
      ignore2 = f2.readline().strip()
      score2 = f2.readline().strip()

      umi = seq2[0:8]
      header1parts[0] = header1parts[0] + "_" + umi
      header2parts[0] = header2parts[0] + "_" + umi

      header1 = " ".join(header1parts)
      header2 = " ".join(header2parts)

      seq2 = seq2[14:]
      score2 = score2[14:]

      readCount = readCount + 1
      if readCount % 100000 == 0:
        logger.info("%d reads processed" % readCount)

      fout1.write(header1 + "\n")
      fout1.write(seq1 + "\n")
      fout1.write(ignore1 + "\n")
      fout1.write(score1 + "\n")

      fout2.write(header2 + "\n")
      fout2.write(seq2 + "\n")
      fout2.write(ignore2 + "\n")
      fout2.write(score2 + "\n")

logger.info("Done!")