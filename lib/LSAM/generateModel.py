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
parser.add_argument('--generate_iteration', action='store_true', help="Generate data by iteration from 1 to 10")
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

lines = []
with open(args.input, 'r') as sr:
  lines = [line for line in sr]
  
outlines = []  
outputDSA = True
outputBorder = True
outputMortality = True
for line in lines:
  if line.startswith(TAG_NAME):
    outlines.append("%s %s\n" % (TAG_NAME, args.name))
  elif line.startswith(TAG_FILE):
    outlines.append("%s %s\n" % (TAG_FILE, args.name))
  elif line.startswith(TAG_START):
    outlines.append("%s %s\n" % (TAG_START, args.startTime))
  elif line.startswith(TAG_FINISH):
    outlines.append("%s %s\n" % (TAG_FINISH, args.endTime))
  elif line.startswith(TAG_DEFDATADEF):
    if args.defdatadef == '-':
      outlines.append(line)
    else:
      outlines.append("%s %s\n" % (TAG_DEFDATADEF, args.defdatadef))
  elif line.startswith(TAG_OPTDATADEF):
    if args.optdatadef == '-':
      outlines.append(line)
    else:
      outlines.append("%s %s\n" % (TAG_OPTDATADEF, args.optdatadef))
  elif line.startswith(TAG_METHODDEF):
    outlines.append("%s %s\n" % (TAG_METHODDEF, args.method))
  elif line.startswith(TAG_TARGET_DSA):
    if args.dsa == '-':
      outlines.append(line)
    else:
      outlines.append("%s %s\n" % (TAG_TARGET_DSA, args.dsa))
      outputDSA = False
  elif line.startswith(TAG_BORDER_STATE):
    if args.border == '-':
      outlines.append(line)
    else:
      outlines.append("%s %s\n" % (TAG_BORDER_STATE, args.border))
      outputBorder = False
  elif line.startswith(TAG_STATE_MORTALITY):
    if args.mortality == '-':
      outlines.append(line)
    else:
      outlines.append("%s %s\n" % (TAG_STATE_MORTALITY, args.mortality))
      outputMortality = False
  elif line.startswith(TAG_OUTPATH):
    if outputDSA and (args.dsa != '-'):
      outlines.append("%s %s\n" % (TAG_TARGET_DSA, args.dsa))
    if outputBorder and (args.border != '-'):
      outlines.append("%s %s\n" % (TAG_BORDER_STATE, args.border))
    if outputMortality and (args.mortality != '-'):
      outlines.append("%s %s\n" % (TAG_STATE_MORTALITY, args.mortality))
    outlines.append(line)
  else:
    outlines.append(line)

if args.generate_iteration:
  for iter in range(1,11):
    outputfile = args.output.replace(".inp", ".%02d.inp" % iter)
    with open(outputfile, 'w') as sw:
      for oline in outlines:
        if oline.startswith(TAG_ORGAN) or oline.startswith(TAG_PATIENT) or oline.startswith(TAG_STATUS):
          oline = oline.replace("Ord1.txt", "Ord%d.txt" % iter)
        sw.write(oline)
else:
  with open(args.output, 'w') as sw:
    for oline in outlines:
      sw.write(oline)
