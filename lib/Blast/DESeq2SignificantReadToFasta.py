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

sequences = {}
for inputfile in inputfiles:
  with open(inputfile, 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='"')
    header = spamreader.next()
    baseMeanIndex = [i for i in range(len(header)) if header[i] == "baseMean"][0] - 1
    for row in spamreader:
      h = {}
      for idx in range(1, baseMeanIndex):
        h[header[idx]] = int(row[idx])
      if row[0] in sequences:
        oldh = sequence[row[0]]
        for k, v in h.items():
          if k in oldh:
            oldh[k] = oldh[k] + h[k]
          else:
            oldh[k] = h[k]
      else:
        sequences[row[0]] = h

with open(outputfile, 'w') as fasta:
  for seq, h in sequences.items():
    count = sum(h.values())
    fasta.write(">%s_%d\n%s\n" %(seq, count, seq));
    
print("%d sequences prepared for blastn \n" % len(sequences))
