import pysam
import argparse
import sys
import logging
import os
import subprocess
from MgfUtils import MgfItem

DEBUG=True
NOT_DEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Shift precursor mass in MGF file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help="Input MGF files")
parser.add_argument('-o', '--output_prefix', action='store', nargs='?', required=NOT_DEBUG, help="Output file prefix")
parser.add_argument('-d', '--shift_dalton', action='store', nargs='?', type=float, default=7.0, required=NOT_DEBUG, help="Shift precursor in Daltons")
parser.add_argument('-p', '--shift_software', action='store', nargs='?', required=NOT_DEBUG, help="Shift precursor software (quality_scope)")
parser.add_argument('-c', '--shift_option_file', action='store', nargs='?', required=NOT_DEBUG, help="Shift precursor option file for quality_scope")

args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/cqs/shengq2/proteomics/20171026_target_decoy_spectra/mgf/Fusion_HCDOT_Human.mgf";
  args.output_prefix = "/scratch/cqs/shengq2/proteomics/20171026_target_decoy_spectra/mgf/Fusion_HCDOT_Human.shift"
  args.shift_dalton = 7.0
  args.shift_software = "/scratch/cqs/shengq2/proteomics/20171026_target_decoy_spectra/software/quality_scope"
  args.shift_option_file = "/scratch/cqs/shengq2/proteomics/20171026_target_decoy_spectra/database/uniprot_human_20170325_gindex.param"

logger = logging.getLogger('shiftPrecursor')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

inputFile = args.output_prefix + ".input"
outputFile = args.output_prefix + ".output"
paramFile = args.output_prefix + ".param"

logger.info("preparing find optimal shift precursor input ...")

with open(args.shift_option_file, 'r') as sr:
  with open(paramFile, 'w') as sw:
    for line in sr:
      if line.startswith("input_path_search"):
        sw.write("input_path_search=%s\n" % inputFile)
      elif line.startswith("output_path_result"):
        sw.write("output_path_result=%s\n" % outputFile)
      else:
        sw.write("%s\n" % line.rstrip())

with open(args.input, 'r') as sr:
  with open(inputFile, 'w') as sw:
      count = 0
      bInScan = False
      for line in sr:
        if line.startswith("BEGIN IONS"):
          bInScan=True
          count = count+1
          if count % 10000 == 0:
            logger.info("%d" % count)
            #break
        elif line.startswith("TITLE="):
          title = line[6:].rstrip()
        elif line.startswith("PEPMASS="):
          precursorMz = float(line[8:].rstrip())
        elif line.startswith("CHARGE=") and bInScan:
          charge = int(line[7])
          precursorMass = (precursorMz - 1.007825035) * charge
        elif line.startswith("END IONS") and bInScan:
          sw.write("%f\t%f\t0.5\t0\t10\t1\n" % (precursorMass, precursorMass + args.shift_dalton))

logger.info("finding optimal shift precursor ...")
proc = subprocess.Popen([args.shift_software, paramFile], cwd=os.path.dirname(args.shift_software), stdout=subprocess.PIPE)
while True:
  line = proc.stdout.readline().rstrip()
  if "The size of query result is" in line:
    proc.kill()
    break
  elif line == "":
    continue
  else:
    print(line.rstrip())

originalFile = args.output_prefix + ".original.mgf"
centerFile = args.output_prefix + ".center.mgf"
optimalFile = args.output_prefix + ".optimal.mgf"

count = 0
with open(outputFile, 'r') as srShift:
  with open(args.input, 'r') as srMgf:
    with open(originalFile, 'w') as swOriginal:
      with open(centerFile, 'w') as swCenter:
        with open(optimalFile, 'w') as swOptimal:
          #write header
          for line in srMgf:
            if(line.startswith("BEGIN IONS")):
              break
            else:
              swOriginal.write("%s\n" % line.rstrip())
              swCenter.write("%s\n" % line.rstrip())
              swOptimal.write("%s\n" % line.rstrip())
              
          mgfItem = MgfItem()
          for line in srShift:
            if line.startswith("Id"):
              continue
            
            parts = line.split(',')
            original = float(parts[1])
            originalCandidates = int(parts[2])
            center = float(parts[3])
            optimal = float(parts[4])
            optimalCandidates = int(parts[5])
            
            mgfItem.read(srMgf)
            
            if originalCandidates == 0:
              continue
          
            title=mgfItem.getTitle()
            charge=mgfItem.getCharge()
            
            mgfItem.write(swOriginal)

            mgfItem.setTitle("CENTER_" + title)
            mgfItem.setPrecursorMz((center + 1.007825035 * charge) / charge)
            mgfItem.write(swCenter)
            
            mgfItem.setTitle("OPTIMAL_" + title)
            mgfItem.setPrecursorMz((optimal + 1.007825035 * charge) / charge)
            mgfItem.write(swOptimal)

            count = count + 1
            if count % 10000 == 0:
              logger.info("Saved %d entries" % count)
              
logger.info("Done! Totally %d entries saved." % count)
