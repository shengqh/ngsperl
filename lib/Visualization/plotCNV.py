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
    #print(locus)
    self.Locus = locus
    parts = locus.split(":")
    #print(parts)
    self.Chromosome = parts[0]
    positions = parts[1].split("-")
    self.Start = int(positions[0])
    self.End = int(positions[1])
    
  def getLocusFileString(self):
    return("%s_%d_%d" % (self.Chromosome, self.Start, self.End))
  
  def str(self):
    return self.Locus
  
class CNVItem(LocusItem):
  def __init__(self, locus, name, sampleCNVMap):
    self.setLocusString(locus)
    self.Name = name
    self.SampleCNVMap = sampleCNVMap
    
  def overlap(self, locus):
    if self.Chromosome != locus.Chromosome:
      return(False)
    if self.End < locus.Start:
      return(False)
    if self.Start > locus.End:
      return(False)
    return(True)
  
  def contains(self, chromosome, position):
    if self.Chromosome != chromosome:
      return(False)
    if self.Start > position:
      return(False)
    if self.End < position:
      return(False)
    return(True)

def readCNVFile(fileName):
  result = []
  samples = []
  with open(fileName, "r") as fin:
    headers = fin.readline().rstrip().split('\t')
    samples = headers[2:]
    for line in fin:
      parts = line.rstrip().split('\t')
      #print(parts)
      sampleCNVMap = {}
      for sampleIndex in range(2, len(parts)):
        cnv = parts[sampleIndex]
        if cnv != "":
          sampleCNVMap[headers[sampleIndex]] = cnv[0:3]
      result.append(CNVItem(parts[0], parts[1], sampleCNVMap))
  return(result, samples)

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Draw bam plot based on peak list.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help="Input bed file")
  parser.add_argument('-b', '--bamListFile', action='store', nargs='?', required=NOT_DEBUG, help="Sample bam file list")
  parser.add_argument('-c', '--cnvFile', action='store', nargs='?', required=NOT_DEBUG, help="CNV file")
  parser.add_argument('-s', '--sizeFactorFile', action='store', nargs='?', required=NOT_DEBUG, help="Sample chromosome size factor file")
  parser.add_argument('-o', '--output', action='store', nargs='?', required=NOT_DEBUG, help="Output folder")
  parser.add_argument('-m', '--minSampleCount', action='store', nargs='?', default=5, help="Minimum sample display in figure")
  parser.add_argument('-n', '--minNormalSampleCount', action='store', nargs='?', default=2, help="Minimum normal sample display in figure")
  parser.add_argument('--geneWithCNVOnly', action='store', nargs='?', default=True, help="Output gene with CNV detected only")

  args = parser.parse_args()
  
  if(DEBUG):
#    args.input = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/20190504_smallRNA_as_chipseq__fileList1.list"
#    args.groupsFile = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/20190504_smallRNA_as_chipseq__fileList2.list"
#    args.bamListFile = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/20190504_smallRNA_as_chipseq__fileList3.list"
#    args.output = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/Control.pdf"
    args.input = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/annotation_genes_locus/result/linton_exomeseq_3321.bed"
    args.bamListFile = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_8_PlotGeneCNV/result/linton_exomeseq_3321__fileList3.list"
    args.cnvFile =  "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_7_CombineGCNV/result/linton_exomeseq_3321.txt"
    args.output = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_8_PlotGeneCNV/result/linton_exomeseq_3321.position.txt"
    args.sizeFactorFile = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/background/linton_exomeseq_3321.excluded.bed.sizefactor"
  
  logger = logging.getLogger('plotPeaks')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  print(args)
  
  bamMap = readHashMap(args.bamListFile)
  cnvs, samples = readCNVFile(args.cnvFile)
  
  #print(samples)

  outputFolder = os.path.dirname(args.output)

  bedFile = args.input
  logger.info("processing " + bedFile + "...")

  bedResultFile = args.output
  bedResultTmpFile = bedResultFile + ".tmp"
  bedResultSlimFile = bedResultFile + ".slim"
  with open(bedResultTmpFile, "wt") as fout:
    fout.write("File\tFeature\tChromosome\tPosition\tPositionCount\tMaxCount\tPercentage\tCNV\n")
    with open(bedResultSlimFile, "wt") as fslim:
      fslim.write("File\tFeature\tCNV\n")

      posData = []
      with open(bedFile, "rt") as fin:
        for line in fin:
          parts = line.rstrip().split('\t')
          chrom = parts[0]
          start = int(parts[1])
          end = int(parts[2])
            
          locus = LocusItem()
          locus.setLocus(chrom, start, end)
            
          overlapCNVs = [cnv for cnv in cnvs if cnv.overlap(locus)]
          
          if (len(overlapCNVs) == 0) and args.geneWithCNVOnly:
            continue

          logger.info("  processing " + locus.Locus + " ...")
  
          if len(parts) > 4:
            locusName = parts[4]
          else:
            locusName = locus.Locus
  
          sampleCNVMap = {}
          for cnv in overlapCNVs:
            for sample in cnv.SampleCNVMap.keys():
              if not sample in sampleCNVMap.keys():
                sampleCNVMap[sample] = [cnv.SampleCNVMap[sample]]
              else:
                sampleCNVMap[sample].append(cnv.SampleCNVMap[sample])
          
          for sample in sorted(sampleCNVMap.keys()):
            sampleCNVs = sorted(set(sampleCNVMap[sample]))
            for cnv in sampleCNVs:
              fslim.write("%s\t%s\t%s\n" % (sample, locusName, cnv))
            
          cnvSamples = set(sampleCNVMap.keys())
          
          otherSamples = sorted(set(samples) - cnvSamples)
          
          if args.minSampleCount > len(cnvSamples):
            normalSampleCount = max(args.minNormalSampleCount, args.minSampleCount - len(cnvSamples))
          else: 
            normalSampleCount = args.minNormalSampleCount
  
          normalSampleCount = min(normalSampleCount, len(otherSamples))
          sampleNames = sorted(cnvSamples.union(otherSamples[0:normalSampleCount]))
          
          sampleFiles = [bamMap[sampleName][0] for sampleName in sampleNames]
      
          bamListFile = bedResultFile + "_" + locus.getLocusFileString() + ".bam.list"
          with open(bamListFile, "w") as flist:
            for sampleFile in sampleFiles:
              flist.write(sampleFile + "\n")
          
          locusData = []
          locusData.append([]) #add chromosome from depth
          locusData.append([]) #add position from depth
          locusData.append([]) #add is position in CNV
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
            inCNV = None
            for cnv in overlapCNVs:
              if cnv.contains(chromosome, position):
                inCNV = cnv
                break
            locusData[2].append(inCNV)
            
            for idx in range(len(sampleNames)):
              locusData[idx+3].append(int(pparts[idx+2]))
          
          chromosomes = locusData[0]
          positions = locusData[1]
          inCNVs = locusData[2]
          for idx in range(len(sampleNames)):
            sampleCount = locusData[idx+3]
            maxCount = max(sampleCount)
  
            if maxCount == 0:
              fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\tNOREAD\n" % (sampleNames[idx], locusName, chromosomes[idx], positions[idx], 0, 0, 0))
              continue
  
            for cIdx in range(len(positions)):
              if sampleCount[cIdx] != 0:
                inCNV = inCNVs[cIdx]
                cnvType = ""
                if inCNV != None:
                  if sampleNames[idx] in inCNV.SampleCNVMap:
                    cnvType = inCNV.SampleCNVMap[sampleNames[idx]]
                fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%s\n" % (sampleNames[idx], locusName, chromosomes[idx], positions[cIdx], sampleCount[cIdx], maxCount, sampleCount[cIdx] * 1.0 / maxCount, cnvType)) 
  
          os.remove(bamListFile)
      
  if os.path.exists(bedResultFile):
    os.remove(bedResultFile)
  os.rename(bedResultTmpFile, bedResultFile)
    
  realpath = os.path.dirname(os.path.realpath(__file__))
  rPath = realpath + "/plotCNV.r"
  cmd = "R --vanilla -f " + rPath + " --args " + bedResultFile + " " + bedResultFile + " " + args.sizeFactorFile
  logger.info(cmd)
  os.system(cmd)

  logger.info("done.")

main()
