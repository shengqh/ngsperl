#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import gzip
from collections import defaultdict
from intervaltree import IntervalTree

def parse_args():
    parser = argparse.ArgumentParser(description='Combine FPKM values from multiple samples with BED file overlaps.')
    parser.add_argument('-i', '--input', required=True, help='Input file list with FPKM file paths and sample names')
    parser.add_argument('-b', '--bed', required=True, help='BED file to find overlaps with (can be gzipped)')
    parser.add_argument('-o', '--output', required=True, help='Output file path for combined FPKM data')
    parser.add_argument('-g', '--genome_coords_cols', default="0,1,2", help='Fallback column indices (0-based) for chr,start,end if named columns not found')
    parser.add_argument('-f', '--fpkm_col', type=int, default=-1, help='Column index (0-based) for FPKM values (default: last column)')
    parser.add_argument('-s', '--skip_lines', type=int, default=0, help='Number of header lines to skip before the header line in FPKM files (default: 0)')
    parser.add_argument('-d', '--delimiter', default='\t', help='Delimiter in FPKM files (default: tab)')
    return parser.parse_args()

def read_file_list(file_path):
    """Read the file list and return a dictionary of sample names to file paths"""
    samples = {}
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split()
                if len(parts) >= 2:
                    file_path = parts[0]
                    sample_name = parts[1]
                    samples[sample_name] = file_path
    return samples

def load_bed_file(bed_file):
    """Load BED file into interval trees by chromosome"""
    trees = {}
    
    # Determine if file is gzipped
    open_func = gzip.open if bed_file.endswith('.gz') else open
    
    with open_func(bed_file, "rt") as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            if line.startswith('V1'):
                continue

            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                
                if chrom not in trees:
                    trees[chrom] = IntervalTree()
                    
                trees[chrom][start:end] = (start, end)
    
    return trees

def find_column_indices(header, delimiter='\t'):
    """Find the indices of the Chr, Start, and End columns in the header"""
    header_parts = header.strip().split(delimiter)
    
    chr_idx, start_idx, end_idx = -1, -1, -1
    
    for i, col in enumerate(header_parts):
        col_lower = col.lower()
        if col_lower == 'chr':
            chr_idx = i
        elif col_lower == 'start':
            start_idx = i
        elif col_lower == 'end':
            end_idx = i
    
    return chr_idx, start_idx, end_idx

def read_fpkm_file(file_path, sample_name, fallback_genome_cols, fpkm_col_index, bed_trees, skip_lines=0, delimiter='\t'):
    """Read a FPKM file, find overlaps with bed file, and return overlapping regions with FPKM values"""
    overlaps = []
    total_peaks = 0
    
    with open(file_path, 'rt') as f:
        # Skip initial lines if specified
        for _ in range(skip_lines):
            next(f)
        
        if file_path.endswith('.bed'):
          chr_idx, start_idx, end_idx = (0,1,2)
          fpkm_idx = -1
        else:
          # Read the header line
          header = next(f)
          
          parts = header.strip().split(delimiter)

          # Try to find column indices from header
          chr_idx, start_idx, end_idx = find_column_indices(header, delimiter)
        
          # Extract FPKM - default to last column if not specified
          fpkm_idx = fpkm_col_index if fpkm_col_index >= 0 else len(parts) - 1

        # Fall back to user-provided indices if not found
        if chr_idx == -1 or start_idx == -1 or end_idx == -1:
            print(f"Warning: Could not find all named columns in {file_path}, falling back to indices: {fallback_genome_cols}")
            chr_idx, start_idx, end_idx = fallback_genome_cols
            
        # Process data lines
        for line in f:
            if not line.startswith('#'):
                parts = line.strip().split(delimiter)
                
                # Extract coordinates
                try:
                    sample_chr = parts[chr_idx]
                    sample_start = int(parts[start_idx])
                    sample_end = int(parts[end_idx])
                except (IndexError, ValueError) as e:
                    print(f"Warning: Error processing line in {file_path}: {e}")
                    continue
                
                total_peaks += 1

                if fpkm_idx == -1:
                    sample_fpkm = -1
                else:
                    sample_fpkm = parts[fpkm_idx]

                # Find overlaps with bed file
                if sample_chr in bed_trees:
                    bed_overlaps = bed_trees[sample_chr][sample_start:sample_end]
                    for overlap in bed_overlaps:
                        bed_start, bed_end = overlap.data
                        # Store all the required information for each overlap
                        overlaps.append((
                            sample_chr,      # bed_chr (same as sample_chr in this case)
                            bed_start,       # bed_start
                            bed_end,         # bed_end
                            sample_name,     # sample_name
                            sample_chr,      # sample_chr
                            sample_start,    # sample_start
                            sample_end,      # sample_end
                            sample_fpkm      # sample_fpkm
                        ))
                
    return overlaps, total_peaks

def main():
    args = parse_args()
    
    # Parse fallback genome coordinate column indices
    fallback_genome_cols = [int(idx) for idx in args.genome_coords_cols.split(',')]
    if len(fallback_genome_cols) != 3:
        raise ValueError("Fallback genome coordinate columns must specify exactly 3 indices for chr, start, end")
    
    # Load BED file into interval trees
    print(f"Loading BED file: {args.bed}")
    bed_trees = load_bed_file(args.bed)
    
    # Read the file list
    samples = read_file_list(args.input)
    if not samples:
        print("No valid samples found in the input file")
        return
    
    # Initialize a list to store all overlaps
    all_overlaps = []
    
    # Dictionary to track total peaks and overlapping peaks per sample
    peak_stats = {}
    
    # Process each FPKM file
    for sample_name, file_path in samples.items():
        print(f"Processing sample: {sample_name}, file: {file_path}")
        
        if not os.path.exists(file_path):
            print(f"Warning: File not found: {file_path}")
            continue
        
        # Read FPKM values from this file and find overlaps with bed file
        overlaps, total_peaks = read_fpkm_file(file_path, sample_name, fallback_genome_cols, args.fpkm_col, bed_trees, 
                                 args.skip_lines, args.delimiter)
        
        # Store total peaks count
        peak_stats[sample_name] = {'total_peaks': total_peaks, 'overlapping_peaks': 0}
        
        # Add to the combined list
        all_overlaps.extend(overlaps)
    
    # Create DataFrame with the new format
    df = pd.DataFrame(all_overlaps, columns=[
        'bed_chr', 'bed_start', 'bed_end', 
        'sample_name', 'sample_chr', 'sample_start', 'sample_end', 'sample_fpkm'
    ])
    
    # Define custom chromosome order
    def chr_order(chr_name):
      chr_name = str(chr_name).lower().replace('chr', '')
      if chr_name.isdigit():
        return int(chr_name)
      elif chr_name == 'x':
        return 23
      elif chr_name == 'y':
        return 24
      elif chr_name == 'm' or chr_name == 'mt':
        return 25
      else:
        return 100  # Place other chromosomes at the end

    # Create sort key column and sort
    df['chr_order'] = df['bed_chr'].apply(chr_order)
    df = df.sort_values(['chr_order', 'bed_start'])
    df = df.drop('chr_order', axis=1)  # Remove the helper column

    # Write to output file
    df.to_csv(args.output, sep='\t', index=False)
    print(f"Overlap data written to: {args.output}")

    # Summarize peak counts per sample
    print("Generating peak count summary...")

    # Get unique peaks for each sample based on sample coordinates
    unique_df = df.drop_duplicates(subset=['sample_name', 'sample_chr', 'sample_start', 'sample_end'])
    
    # Count unique peaks per sample that overlap with BED file
    overlap_counts = unique_df.groupby('sample_name').size().to_dict()
    
    # Update peak_stats with overlapping peak counts
    for sample, count in overlap_counts.items():
        if sample in peak_stats:
            peak_stats[sample]['overlapping_peaks'] = count
    
    # Create summary dataframe
    summary_data = []
    for sample, stats in peak_stats.items():
        summary_data.append({
            'sample_name': sample,
            'total_peaks': stats['total_peaks'],
            'overlapping_peaks': stats['overlapping_peaks'],
            'overlap_percentage': round(100 * stats['overlapping_peaks'] / stats['total_peaks'], 2) if stats['total_peaks'] > 0 else 0
        })
    
    summary_df = pd.DataFrame(summary_data)
    
    # Sort by sample name
    summary_df = summary_df.sort_values('sample_name')
    
    # Write summary to a file in the same directory as the output
    summary_file = os.path.splitext(args.output)[0] + "_peak_counts.txt"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"Peak count summary written to: {summary_file}")

if __name__ == '__main__':
    main()
