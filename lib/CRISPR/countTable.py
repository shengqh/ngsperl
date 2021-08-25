import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
cqsdir = os.path.abspath(os.path.dirname(currentdir) + "/CQS")
sys.path.insert(0,cqsdir) 

import gzip
import logging
import argparse
import subprocess

from FileListUtils import readUniqueHashMap

DEBUG = True
NOT_DEBUG= not DEBUG

if DEBUG:
  sourceFile="/scratch/cqs/shengq2/test/20210311_CRISPRScreen/mageck_count_table/result/CRISPRScreen__fileList1.list"
  outputFile="/scratch/cqs/shengq2/test/20210311_CRISPRScreen/mageck_count_table/result/CRISPRScreen.count.txt"
else:
  parser = argparse.ArgumentParser(description="Combine CRISPR data to count table.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count list file', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output count table file", required=NOT_DEBUG)

  args = parser.parse_args()
  
  print(args)
  
  sourceFile = args.input
  outputFile = args.output

logger = logging.getLogger('countTable')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

countFiles = readUniqueHashMap(sourceFile)

def read(countFile):
  result = {}
  with open(countFile, "rt") as fin:
    fin.readline()
    for line in fin:
      parts = line.rstrip().split('\t')
      result[parts[0]] = parts
  return(result)

counts={}
sgRNAs={}
for sampleName in countFiles.keys():
  countFile = countFiles[sampleName]
  curCounts = read(countFile)
  counts[sampleName]=curCounts
  for sgRNA in curCounts.keys():
    sgRNAs[sgRNA] = curCounts[sgRNA][1]

sampleNames=sorted(countFiles.keys())
sgRNAIndexNameMap = {}
for sgRNAname in sgRNAs.keys():
  sgRNAIndex = int(sgRNAname[6:])
  sgRNAIndexNameMap[sgRNAIndex] = sgRNAname

sgRNAIndecies = sorted(sgRNAIndexNameMap.keys())

with open(outputFile, "wt") as fout:
  fout.write("sgRNA\tGene\t%s\n" % "\t".join(sampleNames))
  for sgRNAindex in sgRNAIndecies:
    sgRNA=sgRNAIndexNameMap[sgRNAindex]
    fout.write("%s\t%s" % (sgRNA, sgRNAs[sgRNA]))
    for sampleName in sampleNames:
      curCounts = counts[sampleName]
      if sgRNA in curCounts:
        fout.write("\t%s" % curCounts[sgRNA][2])
      else:
        fout.write("\t0")
    fout.write("\n")

logger.info("done.")
