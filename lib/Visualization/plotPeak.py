import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
cqsdir = os.path.abspath(os.path.dirname(currentdir) + "/CQS")
sys.path.insert(0,cqsdir) 

import logging
import argparse
import string
import subprocess

from FileListUtils import readHashMap, readUniqueHashMap
from LocusItem import LocusItem, readBedFile

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Draw bam plot based on peak list.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help="Input peak file list")
  parser.add_argument('-g', '--groupsFile', action='store', nargs='?', required=NOT_DEBUG, help="Sample group file")
  parser.add_argument('-b', '--bamListFile', action='store', nargs='?', required=NOT_DEBUG, help="Sample bam file list")
  parser.add_argument('-o', '--output', action='store', nargs='?', required=NOT_DEBUG, help="Output folder")

  args = parser.parse_args()
  
  if(DEBUG):
    args.input = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/20190504_smallRNA_as_chipseq__fileList1.list"
    args.groupsFile = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/20190504_smallRNA_as_chipseq__fileList2.list"
    args.bamListFile = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/20190504_smallRNA_as_chipseq__fileList3.list"
    args.output = "/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result"
  
  logger = logging.getLogger('plotPeaks')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  print(args)
  
  bedMap = readUniqueHashMap(args.input)
  bamMap = readUniqueHashMap(args.bamListFile)
  groupMap = readHashMap(args.groupsFile)

  for bed in bedMap.keys():
    bedFile = bedMap[bed]
    logger.info("processing " + bedFile + "...")

    locusList = readBedFile(bedFile)

    sampleNames = sorted(groupMap[bed])
    sampleFiles = [bamMap[sampleName] for sampleName in sampleNames]

    bamListFile = args.output + "/" + bed + ".bam.list"
    with open(bamListFile, "w") as fout:
      for sampleFile in sampleFiles:
        fout.write(sampleFile + "\n")

    posData = []
    for locus in locusList:
      locusData = []
      locusData.append([]) #add chrom
      locusData.append([]) #add position
      for sampleName in sampleNames:
        locusData.append([])

      logger.info("  processing " + locus.Locus + " ...")
      posData.append([locus, locusData])

      proc = subprocess.Popen(["samtools", "depth", "-f", bamListFile, "-r", locus, "-d", "0"], stdout=subprocess.PIPE)
      for pline in proc.stdout:
        pparts = pline.rstrip().split("\t")
        locusData[0].append(pparts[0])
        locusData[1].append(int(pparts[1]))
        for idx in range(len(sampleNames)):
          locusData[idx+2].append(int(pparts[idx+2]))
 
    bedResultFile = args.output + "/" + bed + ".position.txt"
    with open(bedResultFile, "w") as fout:
      fout.write("File\tFeature\tStrand\tMaxCount\tPositionCount\tPosition\tPercentage\n")
      for pd in posData:
        locus = pd[0]
        locusData = pd[1]
        positions = locusData[1]
        for idx in range(len(sampleNames)):
          sampleCount = locusData[idx+2]
          maxCount = max(sampleCount)

          if maxCount == 0:
            fout.write("%s\t%s\t*\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locus.Name, 0, 0, positions[0], 0))
            continue

          for cIdx in range(len(positions)):
            if sampleCount[cIdx] != 0:
              fout.write("%s\t%s\t*\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locus.Name, maxCount, sampleCount[cIdx], positions[cIdx], sampleCount[cIdx] * 1.0 / maxCount)) 

    os.remove(bamListFile)
    realpath = os.path.dirname(os.path.realpath(__file__))
    rPath = realpath + "/plotPeak.r"
    os.system("R --vanilla -f " + rPath + " --args " + bedResultFile + " " + bedResultFile)


  logger.info("done.")

if __name__ == "__main__":
    main()
