import argparse
import sys
import logging
import os
import csv
import gzip
  
DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Extract barcode and UMI from 10x fastq first read file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-1', '--input1', action='store', nargs='?', help='Input 10x first read file', required=NotDEBUG)
parser.add_argument('-2', '--input2', action='store', nargs='?', help='Input 10x second read file', required=NotDEBUG)
parser.add_argument('-o', '--output1', action='store', nargs='?', help="Output first read file", required=NotDEBUG)
parser.add_argument('-p', '--output2', action='store', nargs='?', help="Output second read file", required=NotDEBUG)
parser.add_argument('-b', '--barcodeFile', action='store', nargs='?', help="Input barcode white list file", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input1 = "/data/cqs/jonathan_brown_data/3804/FASTQ/3804-LD-2_S90_L001_R1_001.fastq.gz"
  args.input2 = "/data/cqs/jonathan_brown_data/3804/FASTQ/3804-LD-2_S90_L001_R2_001.fastq.gz"
  args.output1 = "/data/cqs/jonathan_brown_data/3804/FASTQ/3804-LD-2_S90_L001_R1_001.processed.fastq.gz"
  args.output2 = "/data/cqs/jonathan_brown_data/3804/FASTQ/3804-LD-2_S90_L001_R2_001.processed.fastq.gz"
  args.barcodeFile = "/data/cqs/jonathan_brown_data/3804/Count/3804-LD-2/filtered_feature_bc_matrix/barcodes.tsv.gz"

logger = logging.getLogger('10xFastq')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

barcodes = set()
with gzip.open(args.barcodeFile, 'rt') as fin:
  for barcode in fin:
    barcode = barcode.rstrip()
    barcode = barcode.replace("-1", "")
    barcodes.add(barcode)

logger.info("Total %d barcode in whitelist" % len(barcodes))

tmpFile1 = args.output1 + ".tmp.gz"
tmpFile2 = args.output2 + ".tmp.gz"

fin1count = 0
logger.info("Processing reads ...")
with gzip.open(args.input1, 'rt') as fin1:
  with gzip.open(args.input2, 'rt') as fin2:
    count = 0
    with gzip.open(tmpFile1, "wt") as fout1:
      with gzip.open(tmpFile2, "wt") as fout2:
        while True:
          query = fin1.readline()
          
          if query == "":
            break

          seq = fin1.readline()
          sig = fin1.readline()
          score = fin1.readline()

          fin1count += 1
          if fin1count % 100000 == 0:
            logger.info("processed %d reads ..." % fin1count)

          q2 = fin2.readline()
          seq2 = fin2.readline()
          sig2 = fin2.readline()
          score2 = fin2.readline()

          barcode = seq[:16]
          if not (barcode in barcodes):
            continue

          count = count + 1

          umi = seq[16:26]
          seq = seq[26:]
          score = score[26:]
          query = "@q%d:%s:%s\n" % (count, barcode, umi)

          fout1.write(query)
          fout1.write(seq)
          fout1.write(sig)
          fout1.write(score)

          fout2.write(query)
          fout2.write(seq2)
          fout2.write(sig2)
          fout2.write(score2)
#
#          if count == 1000:
#            break

if os.path.isfile(args.output1):
  os.remove(args.output1)
if os.path.isfile(args.output2):
  os.remove(args.output2)
os.rename(tmpFile1, args.output1)
os.rename(tmpFile2, args.output2)

logger.info("done.")

