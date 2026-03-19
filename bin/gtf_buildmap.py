#!/usr/bin/env python3

import sys
import argparse

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

def merge_intervals(intervals):
  """Merges overlapping intervals and returns the total covered bases."""
  if not intervals:
    return 0
  # Sort intervals by start position
  intervals.sort(key=lambda x: x[0])
  
  merged = []
  curr_start, curr_end = intervals[0]
  
  for next_start, next_end in intervals[1:]:
    if next_start <= curr_end:  # Overlap found
      curr_end = max(curr_end, next_end)
    else:
      merged.append((curr_start, curr_end))
      curr_start, curr_end = next_start, next_end
  merged.append((curr_start, curr_end))
  
  # Calculate total length: (end - start + 1) for each merged segment
  return sum(e - s + 1 for s, e in merged)

def process_gtf(gtf_file, output_file, gene_key='gene_name'):
  """Process GTF file and write map file"""
  genes = {}
  gene_exons = {}
  
  with open(gtf_file, 'r') as f:
    for line in f:
      if line.startswith('#'):
        continue
      
      fields = line.strip().split('\t')
      if len(fields) < 9:
        continue
      
      feature_type = fields[2]
      if feature_type == 'gene':
        chr_name = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        attributes = parse_gtf_attributes(fields[8])
        
        gene_id = attributes['gene_id']
        
        gene_name = attributes.get(gene_key, gene_id)

        if 'gene_type' in attributes:
          gene_biotype = attributes['gene_type']
        elif 'gene_biotype' in attributes:
          gene_biotype = attributes['gene_biotype']
        else:
          print(attributes)
          gene_biotype = ""
        
        length = end - start + 1
        
        genes[gene_id] = {
          'gene_name': gene_name,
          'length': length,
          'chr': chr_name,
          'start': start,
          'end': end,
          'gene_biotype': gene_biotype
        }
      elif feature_type == "exon":
        attributes = parse_gtf_attributes(fields[8])
        gene_id = attributes['gene_id']
        # if not 'exon_number' in attributes:
        #   raise Exception("No exon_number in " + line)
        start = int(fields[3])
        end = int(fields[4])
        if gene_id not in gene_exons:
          gene_exons[gene_id] = []
        gene_exons[gene_id].append((start, end))
  
  # Write output
  with open(output_file, 'w') as f, open(output_file + ".bed", 'w') as bed_f:
    f.write("gene_id\tgene_name\tlength\tchr\tstart\tend\tgene_biotype\n")
    for gene_id, info in sorted(genes.items()):
      if gene_id in gene_exons:
        gene_length = merge_intervals(gene_exons[gene_id])
      else:
        gene_length = info['length']
      f.write(f"{gene_id}\t{info['gene_name']}\t{gene_length}\t{info['chr']}\t{info['start']}\t{info['end']}\t{info['gene_biotype']}\n")
      bed_f.write(f"{info['chr']}\t{info['start']-1}\t{info['end']}\t{gene_id + '_' + info['gene_name']}\n")

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