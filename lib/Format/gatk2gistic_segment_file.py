import argparse
import gzip
import math
import logging

parser = argparse.ArgumentParser(description="Convert GATK segment CNV to gistic2 segmentation file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input list file', required=True)
parser.add_argument('-n', '--no_chr', action='store_true', default=False, help='Remove chr from chromosome')
parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output file", required=True)

args = parser.parse_args()

def convert(name1, file1, no_chr, fout):
  with gzip.open(file1, "rt") as fin:
    for line in fin:
      if line.startswith('#'):
        continue
      
      parts = line.rstrip().split('\t')
      # if parts[4] == '.':
      #   continue
      
      cn = int(parts[9].split(':')[1])
      # if cn == 2:
      #   continue

      start = int(parts[1])
      end = int(parts[7].replace('END=', ''))
      n_markers = math.ceil((end  + 1 - start) / 1000)
      
      if cn == 0:
        segcn = -5
      else:
        segcn = math.log2(cn) - 1
      
      if no_chr:
        chr = parts[0].replace('chr','')
      else:
        chr = parts[0]
      
      segcnstr = "{:.2f}".format(segcn)
      fout.write(f"{name1}\t{chr}\t{start}\t{end}\t{n_markers}\t{segcnstr}\n")
  
logger = logging.getLogger('convert')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

with open(args.output, "wt") as fout:
  with open(args.input, "rt") as fin:
    for line in fin:
      parts = line.rstrip().split('\t')
      logger.info(f"converting {parts[1]} : {parts[0]}")
      convert(parts[1], parts[0], args.no_chr, fout)

logger.info('done')
