import pysam
import sys
import gzip
import os
import logging
import argparse
import xml.etree.ElementTree as ET
import subprocess
from CountXmlUtils import readCountXmlQueryLocationInFeatures

DEBUG = False
NOT_DEBUG= not DEBUG

if DEBUG:
  genomeListFile="/scratch/cqs/shengq2/vickers/20191112_smallRNA_3018-KCV_76_mouse_v4_tRNA_byTiger/data_visualization/bacteria_count/result/HDL_76__fileList1.list"
  databaseListFile = "/scratch/cqs/shengq2/vickers/20191112_smallRNA_3018-KCV_76_mouse_v4_tRNA_byTiger/data_visualization/bacteria_count/result/HDL_76__fileList2.list"
  taskReadFile = "/scratch/cqs/shengq2/vickers/20191112_smallRNA_3018-KCV_76_mouse_v4_tRNA_byTiger/data_visualization/reads_in_tasks/result/HDL_76.NonParallel.TaskReads.csv"
  outputFile="/scratch/cqs/shengq2/vickers/20191112_smallRNA_3018-KCV_76_mouse_v4_tRNA_byTiger/data_visualization/bacteria_count/result/HDL_76.bacteria.count"
else:
  parser = argparse.ArgumentParser(description="Generate smallRNA count from mapped xml.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-g', '--genomeListFile', action='store', nargs='?', help='Input bacteria genome count xml list file', required=NOT_DEBUG)
  parser.add_argument('-d', '--databaseListFile', action='store', nargs='?', help="Original rRNA database count xml list file", required=NOT_DEBUG)
  parser.add_argument('-t', '--taskReadFile', action='store', nargs='?', help="Task read count file", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output count file", required=NOT_DEBUG)

  args = parser.parse_args()
  
  print(args)
  
  genomeListFile = args.genomeListFile
  databaseListFile = args.databaseListFile
  taskReadFile = args.taskReadFile
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

samples = sorted(set(sample for sample in genomeFileMap.keys()))

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
      query_seq = query.get("seq")
      result.setdefault(query_seq, {})[name] = query_count
    
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
        query_seq = query.get("seq")
        result.setdefault(query_seq, {})[name] = query_count

seq_count = [ [seq, sum(result[seq].values())] for seq in result.keys()]

def sortSecond(val): 
    return val[1]

seq_count.sort(key=sortSecond, reverse=True)

with open(outputFile, "wt") as fout:
  fout.write("Sequence\t%s\n" % "\t".join(samples) )
  for query in seq_count:
    query_seq = query[0]
    fout.write(query_seq)
    count_map = result[query_seq]
    for sample in samples:
      fout.write("\t%d" % (count_map[sample] if sample in count_map.keys() else 0))
    fout.write("\n")

summaryFile = outputFile + ".summary"
with open(summaryFile, "wt") as fout:
  fout.write("Sample\tCount\n")
  for sample in samples:
    sample_count = sum(result[seq][sample] if sample in result[seq].keys() else 0 for seq in result.keys())
    fout.write("%s\t%d\n" % (sample, sample_count))

rscript = os.path.realpath(__file__) + ".R"
subprocess.call("R --vanilla -f " + rscript + " --args \"" + summaryFile + "\" \"" + taskReadFile + "\" \"" + summaryFile + "\"", shell=True)

logger.info("done.")
