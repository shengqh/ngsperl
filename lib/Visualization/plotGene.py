import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
cqsdir = os.path.abspath(os.path.dirname(currentdir) + "/CQS")
sys.path.insert(0,cqsdir) 

import logging
import argparse
import string
import subprocess
from LocusItem import LocusItem, readBedFile
from FileListUtils import readUniqueHashMap

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Draw bam plot based on peak list.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help="Input bed file")
  parser.add_argument('-b', '--bamListFile', action='store', nargs='?', required=NOT_DEBUG, help="Sample bam file list")
  parser.add_argument('-s', '--sizeFactorFile', action='store', nargs='?', required=NOT_DEBUG, help="Sample chromosome size factor file")
  parser.add_argument('-e', '--extend_bases', action='store', type=int, default=0, nargs='?', help="Extending X bases before and after coordinates")
  parser.add_argument('-g', '--plot_gene', action='store_true', help="Plot hg38 gene track")
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
  
  logger = logging.getLogger('plotGene')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  print(args)
  
  bamMap = readUniqueHashMap(args.bamListFile)
  sampleNames = sorted(bamMap.keys())
  sampleFiles = [bamMap[sampleName] for sampleName in sampleNames]

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
    fout.write("File\tFeature\tLocus\tPosition\tPositionCount\tMaxCount\tPercentage\n")

    posData = []

    locusList = readBedFile(bedFile)
    for locus in locusList:
      locus.Chromosome = chrMap[locus.Chromosome]
      locusName = locus.getName()
      locusString = locus.getLocusString(args.extend_bases)

      logger.info("  processing " + locus.getLocusString() + " ...")

      locusData = []
      locusData.append([]) #add position from depth
      for sampleName in sampleNames:
        locusData.append([])

      posData.append([locus, locusData])

      proc = subprocess.Popen(["samtools", "depth", "-f", bamListFile, "-r", locusString, "-d", "0"], stdout=subprocess.PIPE)
      for pline in proc.stdout:
        pparts = pline.rstrip().decode("utf-8").split("\t")

        position = int(pparts[1])
        locusData[0].append(position)
          
        for idx in range(len(sampleNames)):
          locusData[idx+1].append(int(pparts[idx+2]))
        
      positions = locusData[0]
      for idx in range(len(sampleNames)):
        sampleCount = locusData[idx+1]
        if len(sampleCount) == 0:
          maxCount = 0
        else:
          maxCount = max(sampleCount)

        if maxCount == 0:
          fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusName, locusString, locus.Start, 0, 0, 0))
          continue

        lastZero = True
        lastPosition = positions[0] - 1
        for cIdx in range(len(positions)):
          curPosition = positions[cIdx]

          if curPosition != lastPosition + 1:
            if not lastZero:
              fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusName, locusString, lastPosition + 1, 0, maxCount, 0))
              lastZero = True

          if sampleCount[cIdx] != 0:
            if lastZero:
              fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusName, locusString, positions[cIdx] - 1, 0, maxCount, 0)) 
            fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusName, locusString, positions[cIdx], sampleCount[cIdx], maxCount, sampleCount[cIdx] * 1.0 / maxCount)) 
            lastZero = False
          else:
            if not lastZero:
              fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusName, locusString, positions[cIdx], 0, maxCount, 0))
            lastZero = True
          lastPosition = curPosition

        fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusName, locusString, positions[len(positions)-1] + 1, 0, maxCount, 0))
  
  if os.path.exists(bedResultFile):
    os.remove(bedResultFile)
    os.remove(bamListFile)
  os.rename(bedResultTmpFile, bedResultFile)
    
  realpath = os.path.dirname(os.path.realpath(__file__))
  #rPath = realpath + "/plotGeneHuman.r" if args.plot_gene else realpath + "/plotGene.r"
  #plotGeneHuman is still under development
  rPath = realpath + "/plotGene.r" if args.plot_gene else realpath + "/plotGene.r"

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

main()
