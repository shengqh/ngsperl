import sys
import os
import logging
import codecs
import argparse
from CountXmlUtils import readCountXmlFeatures

DEBUG = False

if DEBUG:
  inputFile="/scratch/cqs/shengq2/vickers/20180413_smallRNA_1546_RA_3_mouse/host_genome/bowtie1_genome_1mm_NTA_smallRNA_count/result/HFD_Air_PIGR_S20/HFD_Air_PIGR_S20.count.mapped.xml"
  outputFile="/scratch/cqs/shengq2/temp/HFD_Air_PIGR_S20.tRH.count.mapped.xml"
  minLength = 30
  maxLength = 40
else:
  parser = argparse.ArgumentParser(description="Extract tRH reads from mapped xml.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count xml file')
  parser.add_argument('--minLength', action='store', type=int, nargs='?', default = 30, help="Mininum tRH read length", required=NOT_DEBUG)
  parser.add_argument('--maxLength', action='store', type=int, nargs='?', default = 40, help="Maximum tRH read length", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output count xml file")

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  outputFile = args.output
  minLength = args.minLength
  maxLength = args.maxLength

logger = logging.getLogger('filterTRHXml')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def acceptTRNA(featureName):
  return featureName.startswith("tRNA:")
    
logger.info("Reading queries from " + inputFile)
features = readCountXmlFeatures(inputFile, acceptTRNA)

logger.info("%d tRH group read" % len(features))

for fgroup in features:
  fgroup.Queries = [ q for q in fgroup.Queries if len(q.Sequence) >= minLength and len(q.Sequence) <= maxLength ]
  
features = [f for f in features if len(f.Queries) > 0]
logger.info("%d valid tRH group found" % len(features))

queryNames = set()
featureNames = set()
for fgroup in features:
  for q in fgroup.Queries:
    queryNames.add(q.Name)
  for f in fgroup.Features:
    featureNames.add(f.Name)
    
logger.info("%d queries and %d tRHs found" % (len(queryNames), len(featureNames)))

with codecs.open(outputFile, mode="w", encoding="utf-8") as sw:
  with codecs.open(inputFile, mode="r", encoding="utf-8") as sr:
    sw.write(sr.readline())
    
    #output queries
    outputLine = True
    for line in sr:
      if "</queries>" in line:
        sw.write(line)
        break
      
      if "<query" in line:
        qname = line.split("\"")[1]
        outputLine = qname in queryNames
        
      if outputLine:
        sw.write(line)
      
    #output features
    sw.write(sr.readline())
    newGroup = False
    for line in sr:
      if "</subjectResult>" in line:
        sw.write(line)
        break
      
      if "<subjectGroup>" in line:
        newGroup = True
        outputLine = False
        subjectGroupLine = line
        
      if "<subject name" in line:
        fname = line.split("\"")[1]
        outputLine = fname in featureNames
        if outputLine:
          if newGroup:
            sw.write(subjectGroupLine)
            newGroup = False
      
      if outputLine:
        if "<query qname" in line:
          qname = line.split("\"")[1]
          if qname in queryNames:
            sw.write(line)
        else:
          sw.write(line)
        
    #output rest of documents
    for line in sr:
      sw.write(line)
       
#features2 = readCountXmlFeatures(outputFile)

#assert len(features) == len(features2)
