import subprocess
import os.path
import re
import argparse
import csv

parser = argparse.ArgumentParser(description="Transfer DEseq2 significant read to fasta for blastn.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input DEseq2 significant read files, seperated by ","', required=True)
parser.add_argument('-o', '--output', action='store', nargs='?', help='Output annovar _multianno.txt file', required=True)

args = parser.parse_args()

inputfiles=args.input.split(",")
outputfile=args.output

sequences = set()
for inputfile in inputfiles:
  with open(inputfile, 'rb') as csvfile:
    pamreader = csv.reader(csvfile, delimiter=',', quotechar='"')
    header=pamreader.next()
    for row in spamreader:
      sequences.add(row[0])

print("%d sequences prepared for blastn \n" % len(sequences))

print("Transfer DEseq2 significant read to fasta for blastn done.")