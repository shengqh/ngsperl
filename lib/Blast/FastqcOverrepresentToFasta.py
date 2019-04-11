import subprocess
import os.path
import re
import argparse
import csv

parser = argparse.ArgumentParser(description="Transfer FastQC overrepresent read to fasta for blastn.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input FastQC data files, seperated by ","', required=True)
parser.add_argument('-o', '--output', action='store', nargs='?', help='Output fasta file', required=True)

args = parser.parse_args()

inputfiles=args.input.split(",")
outputfile=args.output

sequences = {}
for inputfile in inputfiles:
  with open(inputfile, 'r') as fin:
    for line1 in fin:
      if line1.startswith(">>Overrepresented sequences"):
        break

    for line2 in fin:
      if line2.startswith("#Sequence"):
        continue
      if line2.startswith(">>END_MODULE"):
        break
      parts = line2.split('\t')
      if parts[3].startswith('No Hit'):
        if parts[0] not in sequences:
          sequences[parts[0]] = int(parts[1])
        else:
          sequences[parts[0]] = sequences[parts[0]] + int(parts[1])
      
sorted_sequences = reversed(sorted((value, key) for (key,value) in sequences.items()))
with open(outputfile, "w") as fout:
  for seq in sorted_sequences:
    fout.write(">%s_%d\n%s\n" % (seq[1], seq[0], seq[1]))
