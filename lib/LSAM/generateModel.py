import argparse
import sys
import os

parser = argparse.ArgumentParser(description="Generate LSAM method file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input template file', required=True)
parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output model file")
parser.add_argument('--name', action='store', nargs='?', default="-", help="Method name")
parser.add_argument('--method', action='store', nargs='?', default="-", help="Method file")
parser.add_argument('--dsa', action='store', nargs='?', default="-", help="DSA file")
parser.add_argument('--border', action='store', nargs='?', default="-", help="Border state file")
parser.add_argument('--mortality', action='store', nargs='?', default="-", help="Mortality rate file")
parser.add_argument('--defdatadef', action='store', nargs='?', default="-", help="Default data definition file")
parser.add_argument('--optdatadef', action='store', nargs='?', default="-", help="Optional data definition file")
parser.add_argument('--startTime', action='store', nargs='?', help="Start time")
parser.add_argument('--endTime', action='store', nargs='?', help="End time")

args = parser.parse_args()

TAG_NAME = 'NAME';
TAG_FILE = 'FILE';
TAG_DESCBEG = '<DESC>';
TAG_DESCEND = '</DESC>';
TAG_START = 'START';
TAG_FINISH = 'FINISH';
TAG_DISCARDCNT = 'DISCARDCNT';
TAG_DEFDATADEF = 'DEFDATADEF';
TAG_OPTDATADEF = 'OPTDATADEF';
TAG_LOCMAP = 'LOCMAP';
TAG_GROUPREGION = 'REGIONS';
TAG_ABOMAP = 'ABOMAP';
TAG_BOOSTDEF = 'BOOSTDEF';
TAG_METHODDEF = 'METHODS';
TAG_ACCEPTDEF = 'ACCEPTDEF';
TAG_PARTIALDEF = 'PARTIALDEF';
TAG_SURVIVDEF = 'SURVIVDEF';
TAG_NRDEATHDEF = 'NRDEATHDEF';
TAG_ORGAN = 'ORGAN';
TAG_PATIENT = 'PATIENT';
TAG_STATUS = 'STATUS';
TAG_PAYBACK = 'PAYBACK';
TAG_ANTIGEN = 'ANTIGEN';
TAG_UNANTEQ = 'UNANTEQ';
TAG_OUTPATH = 'OUTPATH';
TAG_EOF = 'EOF';
TAG_TRANSPORT = 'TRANSPORT';
TAG_TARGET_DSA = 'TARGETDSA';
TAG_BORDER_STATE = 'BORDERSTATE';
TAG_STATE_MORTALITY = 'STATEMORTALITY';
TAG_SHARE_PERCENTAGE = 'SHARE_PERCENTAGE';

with open(args.input, 'r') as sr:
  with open(args.output, 'w') as sw:
    outputDSA = True
    outputBorder = True
    outputMortality = True
    for line in sr:
      if line.startswith(TAG_NAME):
        sw.write("%s %s\n" % (TAG_NAME, args.name))
      elif line.startswith(TAG_FILE):
        sw.write("%s %s\n" % (TAG_FILE, args.name))
      elif line.startswith(TAG_START):
        sw.write("%s %s\n" % (TAG_START, args.startTime))
      elif line.startswith(TAG_FINISH):
        sw.write("%s %s\n" % (TAG_FINISH, args.endTime))
      elif line.startswith(TAG_DEFDATADEF):
        if args.defdatadef == '-':
          sw.write(line)
        else:
          sw.write("%s %s\n" % (TAG_DEFDATADEF, args.defdatadef))
      elif line.startswith(TAG_OPTDATADEF):
        if args.optdatadef == '-':
          sw.write(line)
        else:
          sw.write("%s %s\n" % (TAG_OPTDATADEF, args.optdatadef))
      elif line.startswith(TAG_METHODDEF):
        sw.write("%s %s\n" % (TAG_METHODDEF, args.method))
      elif line.startswith(TAG_TARGET_DSA):
        if args.dsa == '-':
          sw.write(line)
        else:
          sw.write("%s %s\n" % (TAG_TARGET_DSA, args.dsa))
          outputDSA = False
      elif line.startswith(TAG_BORDER_STATE):
        if args.border == '-':
          sw.write(line)
        else:
          sw.write("%s %s\n" % (TAG_BORDER_STATE, args.border))
          outputBorder = False
      elif line.startswith(TAG_STATE_MORTALITY):
        if args.mortality == '-':
          sw.write(line)
        else:
          sw.write("%s %s\n" % (TAG_STATE_MORTALITY, args.mortality))
          outputMortality = False
      elif line.startswith(TAG_OUTPATH):
        if outputDSA and (args.dsa != '-'):
          sw.write("%s %s\n" % (TAG_TARGET_DSA, args.dsa))
        if outputBorder and (args.border != '-'):
          sw.write("%s %s\n" % (TAG_BORDER_STATE, args.border))
        if outputMortality and (args.mortality != '-'):
          sw.write("%s %s\n" % (TAG_STATE_MORTALITY, args.mortality))
        sw.write(line)
      else:
        sw.write(line)
 
