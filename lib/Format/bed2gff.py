#Convert bed file to gff file, filtering by name and creating new names based on chr+position.
import sys
import os
import re
import pandas as pd
import argparse

def bed_to_gff(df):
    """Convert BED DataFrame to GFF format DataFrame."""
    gff_df = pd.DataFrame()
    gff_df['seqid'] = df['chr']
    gff_df['source'] = df['name']  # Standard source field, use name
    gff_df['type'] = 'region'  # Default feature type
    gff_df['start'] = df['start'] + 1  # BED is 0-based, GFF is 1-based
    gff_df['end'] = df['end']
    gff_df['score'] = df['score']
    gff_df['strand'] = df['strand']
    gff_df['phase'] = '.'  # Default phase
    return gff_df

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Convert BED file to GFF format with name filtering.')
    parser.add_argument('--bed', '-b', 
                        required=True,
                        help='Input BED file path')
    parser.add_argument('--name', '-n',
                        required=False,
                        help='Name to filter entries')
    parser.add_argument('--output', '-o',
                        required=True,
                        help='Output GFF file path')
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Read and filter BED file
    df = pd.read_csv(args.bed, sep='\t', header=None)
    column_names = ['chr', 'start', 'end', 'name', 'score', "strand"]
    df.columns = column_names + list(df.columns[6:]) if len(df.columns) > 6 else column_names

    if args.name:
      df = df[df['name'] == args.name]
    
    # Create new names
    df['name'] = df['name'] + '_' + df['chr'].astype(str) + '_' + df['start'].astype(str)
    
    # Convert to GFF format
    gff_df = bed_to_gff(df)
    
    # Write GFF file
    gff_df.to_csv(args.output, sep='\t', header=False, index=False)
    print(f"Converted GFF file written to: {args.output}")

if __name__ == '__main__':
    main()
