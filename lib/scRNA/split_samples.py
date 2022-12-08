import pysam
import argparse
import sys
import logging
import os
import pandas as pd

DEBUG = False
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="split bam to subsamples.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input cell classification file', required=NOT_DEBUG)
parser.add_argument('-b', '--bam', action='store', nargs='?', help='Input BAM file', required=NOT_DEBUG)
parser.add_argument('-s', '--sample', action='store', nargs='?', help='Input sample definition file', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output folder", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="/scratch/cqs/alexander_gelbard_projects/20201202_5126_scRNA_split/split_samples/result/COVID/COVID.HTO.csv"
  args.bam="/data/h_vangard_1/alexander_gelbard_data/AG_5126_10X/Count/5126-AG-4/possorted_genome_bam.bam"
  args.sample="/scratch/cqs/alexander_gelbard_projects/20201202_5126_scRNA_split/split_bam/result/fileList_3_COVID.txt"
  args.output="/scratch/cqs/alexander_gelbard_projects/20201202_5126_scRNA_split/split_bam/result"

logger = logging.getLogger('splitSamples')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

cells=pd.read_csv(args.input)
cells=cells.loc[cells['HTO.global'] == 'Singlet']
barcode_dict = dict(zip(cells.iloc[:, 0], cells.HTO))

samples=pd.read_table(args.sample, header=None)
samples_dict = dict(zip(samples.iloc[:, 1], samples.iloc[:, 0]))

sam_dict = {}
index = 0
with pysam.Samfile(args.bam, "rb") as sam:
  header = sam.header
  for read in sam.fetch(until_eof=True):
    index += 1
    if index % 1000000 == 0:
      logger.info("%s : %d" % (read.reference_name, index))

    if read.is_unmapped:
      continue

    if not read.has_tag('CB'):
      continue
    
    barcode = read.get_tag('CB')
    if barcode in barcode_dict:
      hto = barcode_dict[barcode]
      if hto in samples_dict:
        fname = samples_dict[hto]
        if not fname in sam_dict:
          new_file = os.path.join(args.output, fname + ".bam")
          newsam = pysam.AlignmentFile(new_file, "wb", header=header)
          sam_dict[fname] = newsam
        sam_dict[fname].write(read)

for fout in sam_dict.values():
  fout.close()

for fname in samples_dict.values():
  new_file = os.path.join(args.output, fname + ".bam")
  logger.info(f"indexing {new_file}")
  pysam.index(new_file)