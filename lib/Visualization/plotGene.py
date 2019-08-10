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
  def __init__(self):
    self.Locus = ""
    self.Chromosome = ""
    self.Start = 0
    self.End = 0

  def setLocus(self, chromosome, start, end):
    self.Chromosome = chromosome
    self.Start = start
    self.End = end
    self.Locus = "%s:%d-%d" % (chromosome, start, end)
  
  def setLocusString(self, locus):
    self.Locus = locus
    parts = locus.split(":")
    self.Chromosome = parts[0]
    positions = parts[1].split("-")
    self.Start = int(positions[0])
    self.End = int(positions[1])
    
  def getLocusFileString(self):
    return("%s_%d_%d" % (self.Chromosome, self.Start, self.End))
  
  def str(self):
    return self.Locus
  
def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Draw bam plot based on peak list.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help="Input bed file")
  parser.add_argument('-b', '--bamListFile', action='store', nargs='?', required=NOT_DEBUG, help="Sample bam file list")
  parser.add_argument('-s', '--sizeFactorFile', action='store', nargs='?', required=NOT_DEBUG, help="Sample chromosome size factor file")
  parser.add_argument('-o', '--output', action='store', nargs='?', required=NOT_DEBUG, help="Output folder")

  args = parser.parse_args()
  
  if(DEBUG):
#    args.input = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/20190504_smallRNA_as_chipseq__fileList1.list"
#    args.groupsFile = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/20190504_smallRNA_as_chipseq__fileList2.list"
#    args.bamListFile = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/20190504_smallRNA_as_chipseq__fileList3.list"
#    args.output = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/Control.pdf"
    args.input = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/annotation_genes_locus/result/linton_exomeseq_3321.bed"
    args.bamListFile = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_8_PlotGeneCNV/result/linton_exomeseq_3321__fileList3.list"
    args.output = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_8_PlotGeneCNV/result/linton_exomeseq_3321.position.txt"
    args.sizeFactorFile = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/background/linton_exomeseq_3321.excluded.bed.sizefactor"
  
  logger = logging.getLogger('plotPeaks')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  print(args)
  
  bamMap = readHashMap(args.bamListFile)
  sampleNames = sorted(bamMap.keys())
  sampleFiles = [bamMap[sampleName][0] for sampleName in sampleNames]

  outputFolder = os.path.dirname(args.output)

  bedFile = args.input
  logger.info("processing " + bedFile + "...")

  bedResultFile = args.output

  bamListFile = bedResultFile + ".bam.list"
  with open(bamListFile, "w") as flist:
    for sampleFile in sampleFiles:
      flist.write(sampleFile + "\n")

  chrMap = {}
  with open(args.sizeFactorFile, "rt") as fin:
    for line in fin:
      parts = line.rstrip().split('\t')
      chrom = parts[0]
      chromKey = chrom.replace("chr","")
      chrMap[chromKey] = chrom
      chrMap[chrom] = chrom

  #print(chrMap)

  bedResultTmpFile = bedResultFile + ".tmp"
  with open(bedResultTmpFile, "wt") as fout:
    fout.write("File\tFeature\tChromosome\tPosition\tPositionCount\tMaxCount\tPercentage\n")

    posData = []
    with open(bedFile, "rt") as fin:
      for line in fin:
        parts = line.rstrip().split('\t')
        chrom = chrMap[parts[0]]
        start = int(parts[1])
        end = int(parts[2])
            
        locus = LocusItem()
        locus.setLocus(chrom, start, end)
            
        logger.info("  processing " + locus.Locus + " ...")
  
        if len(parts) > 4:
          locusName = parts[4]
        else:
          locusName = locus.Locus
 
        sampleNames = sorted(bamMap.keys())
        sampleFiles = [bamMap[sampleName][0] for sampleName in sampleNames]
      
        locusData = []
        locusData.append([]) #add chromosome from depth
        locusData.append([]) #add position from depth
        for sampleName in sampleNames:
          locusData.append([])
  
        #locusName = locusName + "(%d/%d)" % (len(cnvSamples), len(samples))
        posData.append([locusName, locusData])
  
        proc = subprocess.Popen(["samtools", "depth", "-f", bamListFile, "-r", locus.Locus, "-d", "0"], stdout=subprocess.PIPE)
        for pline in proc.stdout:
          pparts = pline.rstrip().decode("utf-8").split("\t")

          chromosome = pparts[0]
          position = int(pparts[1])
          locusData[0].append(chromosome)
          locusData[1].append(position)
            
          for idx in range(len(sampleNames)):
            locusData[idx+2].append(int(pparts[idx+2]))
          
        chromosomes = locusData[0]
        positions = locusData[1]
        for idx in range(len(sampleNames)):
          sampleCount = locusData[idx+2]
          if len(sampleCount) == 0:
            maxCount = 0
          else:
            maxCount = max(sampleCount)
  
          if maxCount == 0:
            fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusName, chrom, start, 0, 0, 0))
            continue
  
          lastZero = True
          lastPosition = positions[0] - 1
          for cIdx in range(len(positions)):
            curPosition = positions[cIdx]

            if curPosition != lastPosition + 1:
              if not lastZero:
                fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusName, chromosomes[cIdx], lastPosition + 1, 0, maxCount, 0))
                lastZero = True

            if sampleCount[cIdx] != 0:
              if lastZero:
                fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusName, chromosomes[cIdx], positions[cIdx] - 1, 0, maxCount, 0)) 
              fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusName, chromosomes[cIdx], positions[cIdx], sampleCount[cIdx], maxCount, sampleCount[cIdx] * 1.0 / maxCount)) 
              lastZero = False
            else:
              if not lastZero:
                fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusName, chromosomes[cIdx], positions[cIdx], 0, maxCount, 0))
              lastZero = True
            lastPosition = curPosition

          fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusName, chromosomes[0], positions[len(positions)-1] + 1, 0, maxCount, 0))
  
  if os.path.exists(bedResultFile):
    os.remove(bedResultFile)
    os.remove(bamListFile)
  os.rename(bedResultTmpFile, bedResultFile)
    
  realpath = os.path.dirname(os.path.realpath(__file__))
  rPath = realpath + "/plotGene.r"
  cmd = "R --vanilla -f " + rPath + " --args " + bedResultFile + " " + bedResultFile + " " + args.sizeFactorFile
  logger.info(cmd)
  #os.system(cmd)

  logger.info("done.")

main()
