import argparse
import sys
import logging
import os
from asyncore import read

DEBUG=True
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="filter mutect result to keep tumor sample only.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input vcf file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output vcf file", required=NotDEBUG)

args = parser.parse_args()

if DEBUG:
  args.input = "H:/shengquanhu/projects/20190610_Ciombior_ExomeSeq/ID_04-06829.somatic.pass.vcf"
  args.output = "H:/shengquanhu/projects/20190610_Ciombior_ExomeSeq/ID_04-06829.tumor.vcf"

logger = logging.getLogger('filterMutect')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("Reading file %s ..." % args.input)
with open(args.input, "rt") as fin, open(args.output, "wt") as fout:
  tumor_sample_name = ""
  for line in fin:
    if line.startswith("##"):
      if line.startswith("##GATKCommandLine"):
        if not "tumor_sample_name=" in line:
          raise Exception("The file is not mutect format, I cannot find tumor_sample_name in GATKCommandLine: %s" % line);
        tumor_parts = line.split("tumor_sample_name=")[1]
        tumor_sample_name = tumor_parts.split(" ")[0]
      fout.write(line)
    elif line.startswith("#CHROM"):
      if tumor_sample_name == "":
        raise Exception("The file is not mutect format, I cannot find ##GATKCommandLine in %s" % args.input);
      parts = line.rstrip().split("\t")
      format_index = parts.index("FORMAT")
      tumor_index = parts.index(tumor_sample_name)
      logger.info("format_index=%d; tumor_index=%d" % (format_index, tumor_index))
      indecies = [index for index in range(0, len(parts)) if index <= format_index or index == tumor_index]
      logger.info("Indecies=%s" % str(indecies))
      fout.write("%s\n" % "\t".join([parts[index] for index in indecies]))
    else:
      parts = line.rstrip().split("\t")
      fout.write("%s\n" % "\t".join([parts[index] for index in indecies]))
      
logger.info("done.")
    
