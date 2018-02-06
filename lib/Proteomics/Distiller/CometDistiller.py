import pysam
import argparse
import sys
import logging
import os
import subprocess
from networkx.linalg import spectrum

DEBUG=False
NOT_DEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Comet Peptide Spectrum Match Distiller.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help="Input comet xml files")
parser.add_argument('-o', '--output', action='store', nargs='?', required=NOT_DEBUG, help="Output PSM file")

args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/cqs/shengq2/proteomics/20171026_target_decoy_spectra/HCDOT_Human_comet_7_target/result/QExactive_HCDOT_Human.shifted7daltons.center.pep.xml";
  args.output = "/scratch/cqs/shengq2/proteomics/20171026_target_decoy_spectra/HCDOT_Human_comet_7_target/result/QExactive_HCDOT_Human.shifted7daltons.center.pep.xml.tsv"

logger = logging.getLogger('CometPSM')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def getValue(line, key):
  return line.split(key + "=\"",1)[1].split("\"")[0]

logger.info("processing " + args.input)
with open(args.input, 'r') as sr:
  with open(args.output, 'w') as sw:
    sw.write("FileScan\tCharge\tRank\tMatchCount\tXCorr\tIsDecoy\n")
    #sw.write("FileScan\tObservedMass\tCharge\tRank\tPeptide\tCalcMass\tMassDiff\tMatchCount\tXCorr\n")
    count = 0
    bInScan = False
    for line in sr:
      if line.startswith("<spectrum_query"):
        spectrumId = getValue(line, "spectrumNativeID")
        mass = getValue(line, "precursor_neutral_mass")
        charge = getValue(line, "assumed_charge")
        lastrank = 0
      elif line.startswith("<search_hit"):
        rank = getValue(line, "hit_rank")
        peptide = getValue(line, "peptide")
        protein = getValue(line, "protein")
        isDecoy = "REVERSED_" in protein
        calcmass = getValue(line, "calc_neutral_pep_mass")
        massdiff = getValue(line, "massdiff")
        num_matched_peptides = getValue(line, "num_matched_peptides")
      elif line.startswith("<search_score name=\"xcorr\""):
        xcorr = getValue(line, "value")
      elif line.startswith("</search_hit>"):
        if lastrank != rank:
          lastrank = rank
          sw.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (spectrumId, charge, rank, num_matched_peptides, xcorr, isDecoy))
          #sw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (spectrumId, mass, charge, rank, peptide, calcmass, massdiff, num_matched_peptides, xcorr))
          count = count+1
          if count % 10000 == 0:
            logger.info("%d" % count)
      else:
        continue
              
logger.info("Done! Totally %d entries saved." % count)
