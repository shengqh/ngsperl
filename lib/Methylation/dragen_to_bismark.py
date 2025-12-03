#!/usr/bin/env python3
"""
Convert DRAGEN methylation CX report to Bismark coverage format.

This script processes DRAGEN methylation output line-by-line to avoid
memory issues with large files. It filters for CpG context only and
outputs a gzipped Bismark-compatible coverage file.
"""

import gzip
import os
import sys
import logging
import argparse
from collections import defaultdict
from pathlib import Path

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Convert DRAGEN methylation output to Bismark coverage format',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input CX report file (can be gzipped)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output Bismark coverage file (will be gzipped)'
    )
    parser.add_argument(
        '--context',
        default='CG',
        choices=['CG', 'CHG', 'CHH'],
        help='Methylation context to filter (default: CG for CpG)'
    )
    parser.add_argument(
        '--progress-interval',
        type=int,
        default=10000000,
        help='Report progress every N lines'
    )
    return parser.parse_args()


def open_input_file(filepath):
    """Open input file, handling both gzipped and plain text."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    else:
        return open(filepath, 'r')


def process_dragen_to_bismark(input_file, output_file, target_context='CG', progress_interval=10000000):
    """Convert DRAGEN CX report to Bismark coverage format.
    
    Args:
        input_file: Path to input DRAGEN CX report
        output_file: Path to output Bismark coverage file
        target_context: Methylation context to filter (default: CG)
        progress_interval: Report progress every N lines
    
    Returns:
        dict: Statistics about the conversion
    """
    # Statistics
    context_counts = defaultdict(int)
    filtered_count = 0
    line_count = 0
    skipped_lines = 0
    
    logger.info(f"Processing: {input_file}")
    logger.info(f"Output to: {output_file}")
    logger.info(f"Filtering for context: {target_context}")
    
    # Ensure output directory exists
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Process file line by line
    with gzip.open(output_file, 'wt') as out_f:
        with open_input_file(input_file) as in_f:
            for line in in_f:
                line_count += 1
                
                # Skip empty lines
                line = line.strip()
                if not line:
                    skipped_lines += 1
                    continue
                
                # Parse line: chr pos strand meth unmeth context triNuc
                fields = line.split('\t')
                if len(fields) < 6:
                    logger.warning(f"Line {line_count}: insufficient fields ({len(fields)}), skipping")
                    skipped_lines += 1
                    continue
                
                try:
                    chr_name = fields[0]
                    start = fields[1]
                    strand = fields[2]
                    numCs = int(fields[3])
                    numTs = int(fields[4])
                    context = fields[5]
                except (ValueError, IndexError) as e:
                    logger.warning(f"Line {line_count}: parsing error - {e}, skipping")
                    skipped_lines += 1
                    continue
                
                # Count contexts for statistics
                context_counts[context] += 1
                
                # Report progress
                if line_count % progress_interval == 0:
                    logger.info(f"Progress: {chr_name} - {line_count:,} lines, {filtered_count:,} {target_context} sites")
                
                # Filter for target context only
                if context != target_context:
                    continue
                
                filtered_count += 1
                
                # Calculate coverage percentage
                numAll = numCs + numTs
                coverage = (numCs / numAll * 100) if numAll > 0 else 0.0
                end = start
                
                # Write Bismark format: chr start end coverage numCs numTs
                out_f.write(f"{chr_name}\t{start}\t{end}\t{coverage:.6g}\t{numCs}\t{numTs}\n")
    
    # Return statistics
    return {
        'total_lines': line_count,
        'skipped_lines': skipped_lines,
        'context_counts': dict(context_counts),
        'filtered_count': filtered_count,
        'target_context': target_context
    }


def print_statistics(stats, output_file):
    """Print conversion statistics."""
    logger.info("=" * 60)
    logger.info("Conversion completed successfully")
    logger.info("=" * 60)
    logger.info(f"Total lines processed: {stats['total_lines']:,}")
    logger.info(f"Skipped lines: {stats['skipped_lines']:,}")
    logger.info("\nContext distribution:")
    for ctx in sorted(stats['context_counts'].keys()):
        count = stats['context_counts'][ctx]
        percentage = (count / stats['total_lines'] * 100) if stats['total_lines'] > 0 else 0
        logger.info(f"  {ctx:>5s}: {count:>12,} ({percentage:>5.2f}%)")
    logger.info(f"\n{stats['target_context']} sites written: {stats['filtered_count']:,}")
    logger.info(f"Output file: {output_file}")
    
    # Check output file size
    if os.path.exists(output_file):
        size_mb = os.path.getsize(output_file) / (1024 * 1024)
        logger.info(f"Output file size: {size_mb:.2f} MB")



def main():
    """Main entry point."""
    try:
        # Parse arguments
        args = parse_arguments()
        
        # Validate input file exists
        if not os.path.exists(args.input):
            logger.error(f"Input file not found: {args.input}")
            sys.exit(1)
        
        # Ensure output has .gz extension
        if not args.output.endswith('.gz'):
            args.output += '.gz'
            logger.info(f"Output filename adjusted to: {args.output}")
        
        # Process the conversion
        stats = process_dragen_to_bismark(
            args.input,
            args.output,
            args.context,
            args.progress_interval
        )
        
        # Print statistics
        print_statistics(stats, args.output)
        
    except KeyboardInterrupt:
        logger.error("\nProcess interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()

