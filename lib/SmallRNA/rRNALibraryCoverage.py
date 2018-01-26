import sys
import os
import logging
import argparse
import re
import pysam
from DupCountUtils import readDupCountQueries

DEBUG = 1

if DEBUG:
  inputFile="/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/data_visualization/nonhost_library_rRNA_position_vis/result/KCV_3018_77_78_79__fileList1.list"
  countFile="/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/data_visualization/nonhost_library_rRNA_position_vis/result/KCV_3018_77_78_79__fileList2.list"
  outputFile="/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/data_visualization/nonhost_library_rRNA_position_vis/result/KCV_3018_77_78_79_pw.rRNAlib.position"
  species="SILVA_AAQW01000001.5478011.5479546_Pseudomonas_aeruginosa_PACS2,SILVA_AAQW01000001.6267753.6270631_Pseudomonas_aeruginosa_PACS2,SILVA_AB023371.1.1468_Micrococcus_luteus,SILVA_ADCD01000038.9130.12168_Micrococcus_luteus_SK58,SILVA_AARI02000003.303.3117_Listeria_monocytogenes_FSL_F2-515"
else:
  parser = argparse.ArgumentParser(description="Generate tRNA library coverage.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file list')
  parser.add_argument('-c', '--count', action='store', nargs='?', help="Dupcount file list")
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output coverage file")
  parser.add_argument('-s', '--species', action='store', nargs='?', help="Required species id list, sepearte by ','")

  args = parser.parse_args()
  
  print(args)
  
  inputFile=args.input
  outputFile=args.output
  countFile=args.count
  species=args.species

ids = species.split(',')

logger = logging.getLogger('rRNALibraryCoverage')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def readDict(filename):
  result = {}
  with open(filename, "r") as sr:
    for line in sr:
      parts = line.split('\t')
      result[parts[1].strip()] = parts[0]
  return(result)
  
def addCoverage(speciesMap, read, count):
  if not read.reference_name in speciesMap:
    speciesMap[read.reference_name] = {"TotalCount":count, "Coverages":{}}
  else:
    speciesMap[read.reference_name]["TotalCount"] = speciesMap[read.reference_name]["TotalCount"] + count
  indexMap = speciesMap[read.reference_name]["Coverages"]
  for idx in range(read.reference_start, read.reference_end):
    if not idx in indexMap:
      indexMap[idx] = count
    else:
      indexMap[idx] = indexMap[idx] + count

bamfiles = readDict(inputFile)
countfiles = readDict(countFile)

with open(outputFile, "w") as sw:
  sw.write("File\tFeature\tStrand\tTotalCount\tPositionCount\tPosition\tPercentage\n")
  
  for sampleName in sorted(bamfiles.keys()):
    logger.info("processing %s" % (sampleName))
    countfile = countfiles[sampleName]
    dupCountList = readDupCountQueries(countfile, 1, False)
    dupCountMap = {d.Name:d.Count for d in dupCountList}
    
    bamfile = bamfiles[sampleName]
    openmode = "rb" if bamfile.endswith(".bam") else "r"
    
    speciesMap = {}
    with pysam.AlignmentFile(bamfile, openmode) as samfile:
      processed = 0
      for read in samfile.fetch(until_eof=True):
        processed += 1
        if processed % 1000000 == 0:
          logger.info("processed %d" % (processed))
          
        if read.is_unmapped:
          continue;
        
        if read.reference_name in ids:
          addCoverage(speciesMap, read, dupCountMap[read.qname])
    
    for spec in sorted(speciesMap.keys()):
      totalCount = speciesMap[spec]["TotalCount"]
      covs = speciesMap[spec]["Coverages"]
      poslist = sorted(covs.keys())
      for position in poslist:
        sw.write("%s\t%s\t*\t%d\t%d\t%d\t%.2f\n" %(sampleName, spec, totalCount, covs[position], position, covs[position] * 1.0 / totalCount))
  
logger.info("Result has been saved to %s" % outputFile)
