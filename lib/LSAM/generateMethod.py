import argparse
import sys
import os

parser = argparse.ArgumentParser(description="Generate LSAM method file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input template file', required=True)
parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output method file")
parser.add_argument('--originalName', action='store', nargs='?', help="Original variable name")
parser.add_argument('--targetName', action='store', nargs='?', help="Target variable name")

args = parser.parse_args()

with open(args.input, 'r') as sr:
  with open(args.output, 'w') as sw:
    for line in sr:
      newline = line.replace(args.originalName, args.targetName)
      sw.write(newline)
  
