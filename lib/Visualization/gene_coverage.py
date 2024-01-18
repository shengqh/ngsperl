import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
cqsdir = os.path.abspath(os.path.dirname(currentdir) + "/CQS")
sys.path.insert(0,cqsdir) 

import logging
import argparse
import string
import subprocess
import gzip
from LocusItem import LocusItem, readBedFile
from FileListUtils import readUniqueHashMap

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Draw bam plot based on peak list.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-l', '--locus', action='store', nargs='?', required=NOT_DEBUG, help="Input locus, for example chr17:82314868-82317602, or locus file with locus name")
  parser.add_argument('-n', '--name', action='store', nargs='?', required=NOT_DEBUG, help="Input locus name, for example CD7")
  parser.add_argument('-b', '--bam_list_file', action='store', nargs='?', required=NOT_DEBUG, help="Sample bam file list")
  parser.add_argument('--width', action='store', default=3000, type=int, nargs='?', required=NOT_DEBUG, help="Figure width in pixel")
  parser.add_argument('--height', action='store', default=1500, type=int, nargs='?', required=NOT_DEBUG, help="Figure height in pixel")
  parser.add_argument('-o', '--output_folder', action='store', nargs='?', required=NOT_DEBUG, help="Output folder")

  args = parser.parse_args()
  
  if(DEBUG):
    args.locus = "chr17:82314868-82317602"
    args.name = 'CD7'
    args.bam_list_file = "/scratch/cqs/shengq2/ravi_shah_projects/20220107_rnaseq_discovery_hg38/annotation_genes_plot/result/discovery_cohort__fileList1.list"
    args.output_folder = "/scratch/cqs/shengq2/ravi_shah_projects/20220107_rnaseq_discovery_hg38/figures/"
  
  logger = logging.getLogger('plotGene')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  print(args)

  locusName=args.name
  locusString=None
  if os.path.isfile(args.locus):
    with open(args.locus, "rt") as fin:
      for line in fin:
        parts = line.split('\t')
        if parts[4] == locusName:
          locusString = parts[0] + ":" + parts[1] + "-" + parts[2]
          break
    if locusString is None:
      raise Exception(f'No gene {locusName} found in file {args.locus}')
  else:
    locusString = args.locus    

  locusStart=int(locusString.split(':')[1].split('-')[0])
  bamlist=args.bam_list_file
  output_folder=args.output_folder

  bamMap = {}
  with open(bamlist, "rt") as fin:
    for line in fin:
      parts = line.rstrip().split('\t')
      bamMap[parts[1]] = parts[0]
  # bamMap = {
  #   "P_AI120118_V1":"/scratch/cqs/ravi_shah_projects/20220107_rnaseq_discovery_hg38/star_featurecount/result/P_AI120118_V1_Aligned.sortedByCoord.out.bam",
  #   "P_AI120118_V2": "/scratch/cqs/ravi_shah_projects/20220107_rnaseq_discovery_hg38/star_featurecount/result/P_AI120118_V2_Aligned.sortedByCoord.out.bam"}

  logger = logging.getLogger('coverage')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  sampleNames = sorted(bamMap.keys())
  sampleFiles = [bamMap[sampleName] for sampleName in sampleNames]

  bamListFile = output_folder + "/bam.list"
  with open(bamListFile, "w") as flist:
    for sampleFile in sampleFiles:
      flist.write(sampleFile + "\n")

  bedResultFile = output_folder + "/" + locusName + ".depth.txt.gz"
  with gzip.open(bedResultFile, "wt") as fout:
    fout.write("File\tPosition\tPositionCount\tMaxCount\tPercentage\n")
    posData = []

    locusData = []
    locusData.append([]) #add position from depth
    for sampleName in sampleNames:
      locusData.append([])

    posData.append([locusString, locusData])

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
        fout.write("%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], locusStart, 0, 0, 0))
        continue

      lastZero = True
      lastPosition = positions[0] - 1
      for cIdx in range(len(positions)):
        curPosition = positions[cIdx]

        if curPosition != lastPosition + 1:
          if not lastZero:
            fout.write("%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], lastPosition + 1, 0, maxCount, 0))
            lastZero = True

        if sampleCount[cIdx] != 0:
          if lastZero:
            fout.write("%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], positions[cIdx] - 1, 0, maxCount, 0)) 
          fout.write("%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], positions[cIdx], sampleCount[cIdx], maxCount, sampleCount[cIdx] * 1.0 / maxCount)) 
          lastZero = False
        else:
          if not lastZero:
            fout.write("%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], positions[cIdx], 0, maxCount, 0))
          lastZero = True
        lastPosition = curPosition

      fout.write("%s\t%d\t%d\t%d\t%lf\n" % (sampleNames[idx], positions[len(positions)-1] + 1, 0, maxCount, 0))

  realpath = os.path.dirname(os.path.realpath(__file__))
  rPath = realpath + "/gene_coverage.r"

  targetR = bedResultFile + ".r"
  with open(targetR, "wt") as fout:
    fout.write("inputFile<-\"%s\"\n" % bedResultFile)
    fout.write("name<-\"%s\"\n" % locusName)
    fout.write("locus<-\"%s\"\n" % locusString)
    fout.write("width<-%d\n" % args.width)
    fout.write("height<-%d\n\n" % args.height)
    with open(rPath, "r") as fin:
      for line in fin:
        line = line.rstrip()
        fout.write(line + "\n") 

  cmd = "R --vanilla -f " + targetR
  logger.info(cmd)
  os.system(cmd)

  logger.info("done.")

main()