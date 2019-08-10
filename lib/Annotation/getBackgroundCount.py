import sys
import os
import logging
import argparse
import string
import subprocess
import pysam

def readHashMap(fileName):
  result = {}
  with open(fileName, "r") as f:
    for line in f:
      parts = line.rstrip().split('\t')
      if(len(parts) > 1):
        if parts[1] not in result:
          result[parts[1]] = [parts[0]]
        else:
          result[parts[1]].append(parts[0])
  #print(result)
  return(result)

class LocusItem(object):
  def __init__(self, chromosome, start, end, name):
    self.Chromosome = chromosome
    self.Start = start
    self.End = end
    self.Name = name
    self.Locus = "%s:%d-%d" % (chromosome, start, end)
  
def readBedFile(fileName):
  result = []
  with open(fileName, "r") as fin:
    for line in fin:
      parts = line.rstrip().split('\t')
      item = LocusItem(parts[0], int(parts[1]), int(parts[2]), parts[3])
      result.append(item)
  return(result)

def readCNVName(fileName):
  result = []
  with open(fileName, "r") as fin:
    headers = fin.readline().rstrip().split('\t')
    for line in fin:
      parts = line.rstrip().split('\t')
      result.append(parts[1])
  return(result)

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Calculate size factor of sample in each chromosome.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help="Input including range bed file")
  parser.add_argument('-b', '--bamListFile', action='store', nargs='?', required=NOT_DEBUG, help="Sample bam file list")
  parser.add_argument('-e', '--excludeFile', action='store', nargs='?', help="Excluding range bed file")
  parser.add_argument('-o', '--output', action='store', nargs='?', required=NOT_DEBUG, help="Output file")

  args = parser.parse_args()
  
  if(DEBUG):
#    args.input = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/20190504_smallRNA_as_chipseq__fileList1.list"
#    args.groupsFile = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/20190504_smallRNA_as_chipseq__fileList2.list"
#    args.bamListFile = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/20190504_smallRNA_as_chipseq__fileList3.list"
#    args.output = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/Control.pdf"
    args.input = "/scratch/cqs/references/exomeseq/IDT/xgen-exome-research-panel-targetsae255a1532796e2eaa53ff00001c1b3c.slop50.nochr.bed"
    args.bamListFile = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_9_CNVGenesPlot/result/linton_exomeseq_3321__fileList3.list"
    args.excludeFile =  "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_7_CombineGCNV/result/linton_exomeseq_3321.txt"
    args.output = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/background/linton_exomeseq_3321.excluded.bed"
  
  logger = logging.getLogger('getBackgroundCount')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  print(args)

  bamMap = readHashMap(args.bamListFile)
  samples = sorted(bamMap.keys())
  
  if args.input != None:
    bedItems = readBedFile(args.input)
    if args.excludeFile != None:
      cnvNames = set(readCNVName(args.excludeFile))
    else:
      cnvNames = set()
  
    validBedItems = []
    startIndex = 0
    while True:
      while (startIndex < len(bedItems)) and (bedItems[startIndex].Name in cnvNames):
        startIndex = startIndex + 1
    
      if startIndex == len(bedItems):
        break
    
      startChromosome = bedItems[startIndex].Chromosome
    
      endIndex = startIndex + 1
      while (endIndex < len(bedItems)) and (not bedItems[endIndex].Name in cnvNames) and (bedItems[endIndex].Chromosome == startChromosome):
        endIndex = endIndex + 1
    
      if endIndex == len(bedItems):
        validBedItems.append(LocusItem(bedItems[startIndex].Chromosome, bedItems[startIndex].Start, bedItems[-1].End, ""))
        break
    
      validBedItems.append(LocusItem(bedItems[startIndex].Chromosome, bedItems[startIndex].Start, bedItems[endIndex-1].End, ""))
      startIndex = endIndex

    chromosomes = sorted(set(bi.Chromosome for bi in validBedItems))
  else:
    chromosomes = []
  
  with open(args.output, "w") as fout:
    fout.write("Chromosome\tSample\tCount\n")
    for sample in samples:
      bamFile = bamMap[sample][0]

      with pysam.Samfile(bamFile) as samfile:
        logger.info("start counting %s ..." % sample)
        if len(chromosomes) == 0:
          print samfile.references
          curChromosomes = samfile.references
          for chromosome in curChromosomes:
            chromCount = samfile.count(chromosome)
            logger.info("%s ~ %s : %d" % (sample, chromosome, chromCount))
            fout.write("%s\t%s\t%d\n" % (chromosome, sample, chromCount))

        else:
          for chromosome in chromosomes:
            chromBeds = [bi for bi in validBedItems if bi.Chromosome == chromosome]
            chromCount = 0
            for chromBed in chromBeds:
              chromCount = chromCount + samfile.count(chromBed.Chromosome, chromBed.Start, chromBed.End)
            logger.info("%s ~ %s : %d" % (sample, chromosome, chromCount))
            fout.write("%s\t%s\t%d\n" % (chromosome, sample, chromCount))
       
  realpath = os.path.dirname(os.path.realpath(__file__))
  rPath = realpath + "/getBackgroundCount.r"
  cmd = "R --vanilla -f " + rPath + " --args " + args.output + " " + args.output + ".sizefactor"
  logger.info(cmd)
  os.system(cmd)
          
  logger.info("done.")

main()
