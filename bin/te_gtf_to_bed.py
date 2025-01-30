#!/usr/bin/env python3

import argparse
import sys
import logging

def convert(logger, source_gtf, target_bed):
  with open(target_bed, "wt") as fout:
    with open(source_gtf, "rt") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            attr_dict = {key.strip(): value.strip().strip('"') for key, value in 
                         (item.split() for item in attributes.split(";") if item.strip())}
            gene_id = attr_dict.get("gene_id", "")
            family_id = attr_dict.get("family_id", "")
            class_id = attr_dict.get("class_id", "")
            name = f"{gene_id}:{family_id}:{class_id}"
            fout.write(f"{chrom}\t{int(start)-1}\t{end}\t{name}\t{score}\t{strand}\n")
  logger.info(f"Output file: {target_bed}")

def main():
  parser = argparse.ArgumentParser(description="Convert GTF to BED",
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  parser.add_argument('-i', '--input', action='store', help='Input GTF file', required=True)
  parser.add_argument('-o', '--output', action='store', help='Output BED file', required=True)
  
  if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  logger = logging.getLogger('gtf_to_bed')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  convert(logger, source_gtf=args.input, target_bed=args.output)
  
if __name__ == "__main__":
  main()
