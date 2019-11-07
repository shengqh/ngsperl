import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
cqsdir = os.path.abspath(os.path.dirname(currentdir) + "/CQS")
sys.path.insert(0,cqsdir) 

import logging
import argparse
import string
import subprocess
import pysam
from LocusItem import LocusItem, readBedFile
from CNVItem import CNVItem
from FileListUtils import readUniqueHashMap

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
    args.input = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/annotation_genes_locus/result/linton_exomeseq_3321.bed"
    args.bamListFile = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_8_PlotGeneCNV/result/linton_exomeseq_3321__fileList3.list"
    args.cnvFile =  "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_7_CombineGCNV/result/linton_exomeseq_3321.txt"
    args.output = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_8_PlotGeneCNV/result/linton_exomeseq_3321.position.txt"
    args.sizeFactorFile = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/background/linton_exomeseq_3321.excluded.bed.sizefactor"
  
  logger = logging.getLogger('plotCNV')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  print(args)
  
  bamMap = readUniqueHashMap(args.bamListFile)
  cnvs, samples = readCNVFile(args.cnvFile)
  
  outputFolder = os.path.dirname(args.output)

  bedFile = args.input
  logger.info("processing " + bedFile + "...")

  locusList = readBedFile(bedFile)

  bedResultFile = args.output
  bedResultTmpFile = bedResultFile + ".tmp"
  bedResultSlimFile = bedResultFile + ".slim"
  with open(bedResultTmpFile, "wt") as fout:
    fout.write("File\tFeature\tLocus\tPosition\tPositionCount\tMaxCount\tPercentage\tCNV\n")
    with open(bedResultSlimFile, "wt") as fslim:
      fslim.write("File\tFeature\tCNV\n")

      posData = []

      for locus in locusList:
        overlapCNVs = [cnv for cnv in cnvs if cnv.overlap(locus)]

        locusName = locus.getName()
        locusString = locus.getLocusString()
        
        if (len(overlapCNVs) == 0) and args.geneWithCNVOnly:
          continue

        logger.info("  processing " + locusString + " ...")

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
        
        sampleFiles = [bamMap[sampleName] for sampleName in sampleNames]
    
        bamListFile = bedResultFile + "_" + locus.getLocusFileString() + ".bam.list"
        with open(bamListFile, "w") as flist:
          for sampleFile in sampleFiles:
            flist.write(sampleFile + "\n")
        
        locusData = []
        locusData.append([]) #add position from depth
        locusData.append([]) #add if position is in CNV
        for sampleName in sampleNames:
          locusData.append([])

        #locusName = locusName + "(%d/%d)" % (len(cnvSamples), len(samples))
        posData.append([locusName, locusData])

        proc = subprocess.Popen(["samtools", "depth", "-f", bamListFile, "-r", locusString, "-d", "0"], stdout=subprocess.PIPE)
        for pline in proc.stdout:
          pparts = pline.rstrip().decode("utf-8").split("\t")

          chromosome = pparts[0]
          position = int(pparts[1])
          locusData[0].append(position)
          inCNV = None
          for cnv in overlapCNVs:
            if cnv.contains(chromosome, position):
              inCNV = cnv
              break
          locusData[1].append(inCNV)
          
          for idx in range(len(sampleNames)):
            locusData[idx+2].append(int(pparts[idx+2]))
        
        positions = locusData[0]
        inCNVs = locusData[1]
        for idx in range(len(sampleNames)):
          sampleCount = locusData[idx+2]
          maxCount = max(sampleCount)

          if maxCount == 0:
            fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\tNOREAD\n" % (sampleNames[idx], locusName, locusString, positions[idx], 0, 0, 0))
            continue

          lastZero = True
          lastPosition = positions[0] - 1
          for cIdx in range(len(positions)):
            curPosition = positions[cIdx]
            inCNV = inCNVs[cIdx]
            cnvType = ""
            if inCNV != None:
              if sampleNames[idx] in inCNV.SampleCNVMap:
                cnvType = inCNV.SampleCNVMap[sampleNames[idx]]

            if curPosition != lastPosition + 1:
              if not lastZero:
                fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%s\n" % (sampleNames[idx], locusName, locusString, positions[cIdx], sampleCount[cIdx], maxCount, sampleCount[cIdx] * 1.0 / maxCount, cnvType)) 
                fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%s\n" % (sampleNames[idx], locusName, locusString, lastPosition + 1, 0, maxCount, 0, ""))
                lastZero = True

            if sampleCount[cIdx] != 0:
              if lastZero:
                fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%s\n" % (sampleNames[idx], locusName, locusString, positions[cIdx] - 1, 0, maxCount, 0, cnvType))
              fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%s\n" % (sampleNames[idx], locusName, locusString, positions[cIdx], sampleCount[cIdx], maxCount, sampleCount[cIdx] * 1.0 / maxCount, cnvType)) 
              lastZero = False
            else:
              if not lastZero:
                fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%s\n" % (sampleNames[idx], locusName, locusString, positions[cIdx], 0, maxCount, 0, cnvType))
              lastZero = True
            lastPosition = curPosition

          fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%s\n" % (sampleNames[idx], locusName, locusString, positions[len(positions)-1] + 1, 0, maxCount, 0, cnvType))

        os.remove(bamListFile)
      
  if os.path.exists(bedResultFile):
    os.remove(bedResultFile)
  os.rename(bedResultTmpFile, bedResultFile)
    
  realpath = os.path.dirname(os.path.realpath(__file__))
  rPath = realpath + "/plotGene.r"

  targetR = bedResultFile + ".r"
  with open(targetR, "wt") as fout:
    fout.write("inputFile<-\"%s\"\n" % bedResultFile)
    fout.write("outputPrefix<-\"%s\"\n" % bedResultFile)
    fout.write("sizeFactorFile<-\"%s\"\n\n" % args.sizeFactorFile)
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
