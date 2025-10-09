import os
import argparse
import pysam

def change_quality_scores(bam_file, output_file):
  """
  Changes the quality scores of all reads in a BAM file to F.

  Args:
    bam_file: Path to the input BAM file.
    output_file: Path to the output BAM file.
  """

  # Open input and output BAM files
  with pysam.AlignmentFile(bam_file, "rb") as input_samfile:
    # Get header and remove lines containing ID:DNMTOOLS
    header = input_samfile.header.to_dict()
    if 'PG' in header:
      idx = 0
      # Modify DNMTOOLS ID if present
      for pg in header['PG']:
        if 'ID' in pg and pg['ID'] == 'DNMTOOLS':
          idx = idx + 1
          pg['ID'] = 'DNMTOOLS.' + str(idx)
    with pysam.AlignmentFile(output_file, "wb", header=header) as output_samfile:
      # Iterate over reads in the input BAM file
      for read in input_samfile:
        read.query_qualities  = pysam.qualitystring_to_array("F" * read.query_length)
        output_samfile.write(read)

parser = argparse.ArgumentParser(description='add quality columns to bam files')
parser.add_argument('bam_file', type=str, help='input bam file')
parser.add_argument('output_file', type=str, help='output bam file after fixing quality column')

args = parser.parse_args()
# Run the function
change_quality_scores(args.bam_file, args.output_file)

#print(f"Quality scores successfully changed to F in {output_file}")

