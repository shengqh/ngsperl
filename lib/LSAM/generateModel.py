import argparse
import sys
import os

parser = argparse.ArgumentParser(description="Generate LSAM method file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input template file', required=True)
parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output model file")
parser.add_argument('--name', action='store', nargs='?', default="-", help="Method name")
parser.add_argument('--method', action='store', nargs='?', default="-", help="Method file")
parser.add_argument('--datadef', action='store', nargs='?', default="-", help="Optional data definition file")
parser.add_argument('--startTime', action='store', nargs='?', help="Start time")
parser.add_argument('--endName', action='store', nargs='?', help="End time")

args = parser.parse_args()

with open(args.input, 'r') as sr:
  with open(args.output, 'w') as sw:
    for line in sr:
      if line.startswith("NAME"):
        sw.write("NAME " + args.name + "\n")
      elif line.startswith("FILE"):
        sw.write("FILE " + args.name + "\n")
      elif line.startswith("START"):
        sw.write("START " + args.startTime + "\n")
      elif line.startswith("FINISH"):
        sw.write("FINISH " + args.endName + "\n")
      elif line.startswith("OPTDATADEF"):
        sw.write("OPTDATADEF " + args.datadef + "\n")
      elif line.startswith("METHODS"):
        sw.write("METHODS " + args.method + "\n")
      else:
        sw.write(line)
  
