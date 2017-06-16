import pysam
import argparse
import sys
import logging
import os
from asyncore import read
from Bio import SeqIO

def prepare(output, reference, referencePrefix, homology, homologyPrefix):
  with open(output, "w") as out:
    for record in SeqIO.parse(open(reference),'fasta'):
      record.id = referencePrefix + record.id
      record.name = ""
      record.description = ""
      SeqIO.write(record, out, "fasta")
    for record in SeqIO.parse(open(homology),'fasta'):
      record.id = homologyPrefix + record.id
      record.name = ""
      record.description = ""
      SeqIO.write(record, out, "fasta")
        
def main():
  parser = argparse.ArgumentParser(description="Prepare homology smallRNA database.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = True
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-r', '--reference', action='store', nargs='?', help='Input reference FASTA file)', required=NOT_DEBUG)
  parser.add_argument('--referencePrefix', action='store', nargs='?', help='Input reference genome prefix)', required=NOT_DEBUG)
  parser.add_argument('-l', '--homology', action='store', nargs='?', help="Input homology FASTA file", required=NOT_DEBUG)
  parser.add_argument('--homologyPrefix', action='store', nargs='?', help='Input homology genome prefix)', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output FASTA file", required=NOT_DEBUG)
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  
  if DEBUG:
    args.reference="/scratch/cqs/shengq1/references/smallrna/v3/mm10_miRBase21_GtRNAdb2_gencode12_ncbi.fasta"
    args.referencePrefix="mm10_"
    args.homology="/scratch/cqs/shengq1/references/smallrna/v3/rn5_miRBase21_GtRNAdb2_ensembl79_ncbi.fasta"
    args.homologyPrefix="rn5_"
    args.output="/scratch/cqs/shengq1/references/smallrna/v3/homology/mm10_rn5.fa"
  
  logger = logging.getLogger('homology')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  prepare(args.output, args.reference, args.referencePrefix, args.homology, args.homologyPrefix)
  
if __name__ == "__main__":
    main()
