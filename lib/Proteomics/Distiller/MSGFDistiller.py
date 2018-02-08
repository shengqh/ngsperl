import math
import argparse
import sys
import logging
import os
import subprocess

DEBUG=False
NOT_DEBUG=not DEBUG

parser = argparse.ArgumentParser(description="MSGF Peptide Spectrum Match Distiller.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help="Input MSGF xml files")
parser.add_argument('-o', '--output', action='store', nargs='?', required=NOT_DEBUG, help="Output PSM file")

args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/cqs/shengq2/proteomics/20171026_target_decoy_spectra/HCDOT_Human_MSGF_7_decoy/result/QExactive_HCDOT_Human.shifted7daltons.optimal.msgf.mzid";
  args.output = "/scratch/cqs/shengq2/proteomics/20171026_target_decoy_spectra/HCDOT_Human_MSGF_7_decoy/result/QExactive_HCDOT_Human.shifted7daltons.optimal.msgf.mzid.tsv";

logger = logging.getLogger('MSGF_PSM')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def getValue(line, key):
  return line.split(key + "=\"",1)[1].split("\"")[0]

logger.info("processing " + args.input)
with open(args.input, 'r') as sr:
  with open(args.output, 'w') as sw:
    sw.write("FileScan\tCharge\tRank\tScore\n")
    count = 0
    bInScan = False
    for line in sr:
      if "<SpectrumIdentificationResult" in line:
        spectra={}
      elif "<SpectrumIdentificationItem" in line:
        rank = int(getValue(line, "rank"))
        spectrum = {}
        spectra[rank] = spectrum
        spectrum["charge"] = getValue(line, "chargeState")
      elif "MS-GF:SpecEValue" in line:
        spectrum["score"] = float(getValue(line, "value"))
      elif "spectrum title" in line:
        spectrumId = getValue(line, "value")
      elif "</SpectrumIdentificationResult>" in line:
        for rank in sorted(spectra.keys()):
          spectrum = spectra[rank]
          score = spectrum["score"]
          if(rank > 1 and score == spectra[rank-1]["score"]):
            continue
          sw.write("%s\t%s\t%s\t%s\n" % (spectrumId, spectrum["charge"], rank, -math.log(score)))
        count = count+1
        if count % 10000 == 0:
          logger.info("%d" % count)
      else:
        continue
              
logger.info("Done! Totally %d entries saved." % count)
