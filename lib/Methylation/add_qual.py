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
    with pysam.AlignmentFile(output_file, "wb", template=input_samfile) as output_samfile:
      # Iterate over reads in the input BAM file
      for read in input_samfile:
        # Get the read length
        read_length = len(read.query_sequence)

        # Create an array of F characters with the same length as the read
        new_quality_scores = pysam.qualitystring_to_array("F" * read_length)

        # Update the read's quality scores with the new scores
        read.query_qualities  = new_quality_scores

        # Write the modified read to the output BAM file
        output_samfile.write(read)

parser = argparse.ArgumentParser(description='add quality columns to bam files')
parser.add_argument('bam_file', type=str, help='input bam file')
parser.add_argument('output_file', type=str, help='output bam file after fixing quality column')

args = parser.parse_args()
# Run the function
change_quality_scores(args.bam_file, args.output_file)

#print(f"Quality scores successfully changed to F in {output_file}")

