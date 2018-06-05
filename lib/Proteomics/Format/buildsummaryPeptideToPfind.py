import argparse
import sys
import logging
import os
import subprocess
import csv

DEBUG=True
NOT_DEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Convert buildsummary peptides to pFind peptides file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help="Input buildsummary peptides file")
parser.add_argument('-o', '--output', action='store', nargs='?', required=NOT_DEBUG, help="Output pFind peptides file")

args = parser.parse_args()

if DEBUG:
  args.input = "H:/shengquanhu/projects/20170503_haake_proteogenomics/Proteomics/1DLC/MaxQuant_1.5.8.3_novariants/txt/allPeptides.txt";
  args.output = "H:/shengquanhu/projects/20170503_haake_proteogenomics/Proteomics/1DLC/pFind.spectra";

logger = logging.getLogger('converter')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

with open(args.output, 'w') as sw:
  sw.write("File_Name\tScan_No\tExp.MH+\tCharge\tQ-value\tSequence\tCalc.MH+\tMass_Shift(Exp.-Calc.)\tRaw_Score\tFinal_Score\tModification\tSpecificity\tProteins\tPositions\tLabel\tTarge/Decoy\tMiss.Clv.Sites\tAvg.Frag.Mass.Shift\tOthers\n")
  with open(args.input, 'r') as csvfile:
    mycsv = csv.DictReader(csvfile, delimiter="\t")
    for row in mycsv:
      scan = row["MSMS Scan Numbers"].split(';')[0]
      sw.write("%s.%s.%s.%s.0.dta\t%s\tExp.MH+\tCharge\tQ-value\tSequence\tCalc.MH+\tMass_Shift(Exp.-Calc.)\tRaw_Score\tFinal_Score\tModification\tSpecificity\tProteins\tPositions\tLabel\tTarge/Decoy\tMiss.Clv.Sites\tAvg.Frag.Mass.Shift\tOthers\n" %
               (row["Raw file"], scan, scan, row["Charge"],
                scan,
                ))
                    
logger.info("Done! Totally %d entries saved." % count)
