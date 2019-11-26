import pysam
import sys
import gzip
import os
import logging
import argparse
import xml.etree.ElementTree as ET
from CountXmlUtils import readCountXmlQueryLocationInFeatures

DEBUG = True
NOT_DEBUG= not DEBUG

if DEBUG:
  genomeListFile="/scratch/cqs/shengq2/vickers/20191112_smallRNA_3018-KCV_76_mouse_v4_tRNA_byTiger/data_visualization/bacteria_count/result/HDL_76__fileList1.list"
  databaseListFile = "/scratch/cqs/shengq2/vickers/20191112_smallRNA_3018-KCV_76_mouse_v4_tRNA_byTiger/data_visualization/bacteria_count/result/HDL_76__fileList2.list"
  outputFile="/scratch/cqs/shengq2/vickers/20191112_smallRNA_3018-KCV_76_mouse_v4_tRNA_byTiger/data_visualization/bacteria_count/result/HDL_76.bacteria.count"
else:
  parser = argparse.ArgumentParser(description="Generate smallRNA BAM from mapped xml.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-g', '--genomeListFile', action='store', nargs='?', help='Input bacteria genome count xml list file', required=NOT_DEBUG)
  parser.add_argument('-d', '--databaseListFile', action='store', nargs='?', help="Original rRNA database count xml list file", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output bam file", required=NOT_DEBUG)

  args = parser.parse_args()
  
  print(args)
  
  genomeListFile = args.genomeListFile
  databaseListFile = args.databaseListFile
  outputFile = args.output

logger = logging.getLogger('getBacteriaCount')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def readFileMap(fileName):
  result = {}
  with open(fileName) as fh:
    for line in fh:
      filepath, name = line.strip().split('\t', 1)
      result.setdefault(name, []).append(filepath.strip())
  return(result)
   
genomeFileMap = readFileMap(genomeListFile)
dbFileMap = readFileMap(databaseListFile)

result = {}
for name in genomeFileMap.keys():
  logger.info("Parsing " + name)
  queryMap = {}

  for genomeFile in genomeFileMap[name]:
    logger.info("  Parsing " + genomeFile)
    tree = ET.parse(genomeFile)
    root = tree.getroot()
    queries = root.find('queries')
    for query in queries.findall('query'):
      query_count = int(query.get("count"))
      query_name = query.get("name")
      queryMap[query_name] = query_count
    
  for dbFile in dbFileMap[name]:
    logger.info("  Parsing " + dbFile)
    tree = ET.parse(dbFile)
    root = tree.getroot()
    queries = root.find('queries')
    for query in queries.findall('query'):
      is_bacteria = False
      for loc in query.findall('location'):
        seqname = loc.get("seqname")
        if seqname == "Bacteria":
          is_bacteria = True
          break
      
      if is_bacteria:
        query_count = int(query.get("count"))
        query_name = query.get("name")
        queryMap[query_name] = query_count
        
  result[name] = sum(queryMap.values())   

with open(outputFile, "wt") as fout:
  fout.write("Sample\tCount\n")
  for name in sorted(result.keys()):
    fout.write("%s\t%d\n" % (name, result[name]))
    
logger.info("done.")
