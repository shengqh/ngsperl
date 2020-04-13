import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
cqsdir = os.path.abspath(os.path.dirname(currentdir) + "/CQS")
sys.path.insert(0,cqsdir) 

import logging
import argparse
import string
import subprocess
import pysam

from LocusItem import LocusItem, readBedFile, getChromosomeMap
from FileListUtils import readHashMap 

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG

  parser = argparse.ArgumentParser(description="Draw bam plot based on peak list.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help="Input bed file")
  parser.add_argument('-b', '--bamListFile', action='store', nargs='?', required=NOT_DEBUG, help="Sample bam file list")
  parser.add_argument('-c', '--cnvFile', action='store', nargs='?', required=NOT_DEBUG, help="Exclude CNV range file")
  parser.add_argument('-o', '--output', action='store', nargs='?', required=NOT_DEBUG, help="Output file")

  args = parser.parse_args()

  if(DEBUG):
    args.input = "/scratch/cqs/references/exomeseq/IDT/xgen-exome-research-panel-targetsae255a1532796e2eaa53ff00001c1b3c.slop50.nochr.bed"
    args.bamListFile = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_9_CNVGenesPlot/result/linton_exomeseq_3321__fileList3.list"
    args.cnvFile =  "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_7_CombineGCNV/result/linton_exomeseq_3321.txt"
    args.output = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/background/linton_exomeseq_3321.excluded.bed"

  logger = logging.getLogger('getBackgroundCount')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  print(args)

  #if not os.path.isfile(args.output):
  bamMap = readHashMap(args.bamListFile)
  samples = sorted(bamMap.keys())

  bedItems = readBedFile(args.input)
  cnvItems = readBedFile(args.cnvFile)

  bedMap = getChromosomeMap(bedItems)
  cnvMap = getChromosomeMap(cnvItems)

  logger.info("Before excluding, there are %d intervals" % len(bedItems))
  for chrom in bedMap.keys():
    logger.info(chrom)
    if not chrom in cnvMap:
      continue
    curBedItems = bedMap[chrom]
    curExcludeItems = cnvMap[chrom]
      
    for bi in curBedItems:
      for ci in curExcludeItems:
        if bi.overlapPosition(ci):
          bi.Overlapped = True
          break

  bedItems = [bi for bi in bedItems if not bi.Overlapped]
  logger.info("After excluding, there are %d intervals" % len(bedItems))

  for chrom in bedMap.keys():
    curBedItems = bedMap[chrom]
      
    for idx in range(len(curBedItems)-1, 1, -1):
      curItem = curBedItems[idx]
      prevItem = curBedItems[idx - 1]
      if (not curItem.Overlapped) and (not prevItem.Overlapped):
        prevItem.End = curItem.End
        curItem.Overlapped = True
  
  validBedItems = [bi for bi in bedItems if not bi.Overlapped]
  logger.info("After merge, there are %d intervals" % len(validBedItems))

  chromosomes = sorted(set(bi.Chromosome for bi in validBedItems))
  print(chromosomes)

  with open(args.output, "w") as fout:
    fout.write("Chromosome\tSample\tCount\n")
    for sample in samples:
      bamFile = bamMap[sample][0]

      with pysam.Samfile(bamFile) as samfile:
        logger.info("start counting %s ..." % sample)
        for chromosome in chromosomes:
          chromBeds = [bi for bi in validBedItems if bi.Chromosome == chromosome]
          chromCount = 0
          for chromBed in chromBeds:
            chromCount = chromCount + samfile.count(chromBed.Chromosome, chromBed.Start, chromBed.End)
          logger.info("%s ~ %s : %d" % (sample, chromosome, chromCount))
          fout.write("%s\t%s\t%d\n" % (chromosome, sample, chromCount))

  realpath = os.path.dirname(os.path.realpath(__file__))
  rPath = realpath + "/getBackgroundCount.r"

  targetR = args.output + ".r"
  with open(targetR, "wt") as fout:
    fout.write("inputFile<-\"%s\"\n" % args.output)
    fout.write("outputFile<-\"%s\"\n\n" % (args.output + ".sizefactor"))
    with open(rPath, "r") as fin:
      for line in fin:
        line = line.rstrip()
        fout.write(line + "\n") 

  cmd = "R --vanilla -f " + targetR
  logger.info(cmd)
  os.system(cmd)

  logger.info("done.")

if __name__ == "__main__":
    main()
