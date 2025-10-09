#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict

def parse_gtf_attributes(attr_string):
  """Parse GTF attributes string and return a dictionary"""
  attributes = {}
  for item in attr_string.strip().split(';'):
    if item.strip():
      key_value = item.strip().split(' ', 1)
      if len(key_value) == 2:
        key = key_value[0]
        value = key_value[1].strip('"')
        attributes[key] = value
  return attributes

def process_gtf(gtf_file, output_file, gene_key='gene_name'):
  """Process GTF file and write map file"""
  genes = {}
  
  with open(gtf_file, 'r') as f:
    for line in f:
      if line.startswith('#'):
        continue
      
      fields = line.strip().split('\t')
      if len(fields) < 9:
        continue
      
      feature_type = fields[2]
      if feature_type != 'gene':
        continue
      
      chr_name = fields[0]
      start = int(fields[3])
      end = int(fields[4])
      attributes = parse_gtf_attributes(fields[8])
      
      gene_id = attributes['gene_id']
      
      gene_name = attributes.get(gene_key, gene_id)

      gene_biotype = attributes['gene_biotype']
      
      length = end - start + 1
      
      genes[gene_id] = {
        'gene_name': gene_name,
        'length': length,
        'chr': chr_name,
        'start': start,
        'end': end,
        'gene_biotype': gene_biotype
      }
  
  # Write output
  with open(output_file, 'w') as f:
    f.write("gene_id\tgene_name\tlength\tchr\tstart\tend\tgene_biotype\n")
    for gene_id, info in genes.items():
      f.write(f"{gene_id}\t{info['gene_name']}\t{info['length']}\t{info['chr']}\t{info['start']}\t{info['end']}\t{info['gene_biotype']}\n")

def main():
  parser = argparse.ArgumentParser(description='Build gene map from GTF file')
  parser.add_argument('-i', '--input', help='Input GTF file')
  parser.add_argument('-o', '--output', help='Output map file')
  parser.add_argument('-k', '--gene_key', default='gene_name', help='Key for gene name in GTF attributes')
  args = parser.parse_args()

  process_gtf(args.input, args.output, gene_key=args.gene_key)
  print(f"Gene map written to {args.output}")

if __name__ == '__main__':
  main()