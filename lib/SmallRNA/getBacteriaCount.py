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
  genomeListFile="/scratch/stein_lab/shengq2/20200226_4233_4263_michelle_smallRNA_human_v5_byTiger/data_visualization/bacteria_count/result/StaRRA_human_4233_4263__fileList1.list"
  databaseFile = "/scratch/stein_lab/shengq2/20200226_4233_4263_michelle_smallRNA_human_v5_byTiger/nonhost_library/bowtie1_rRNA_pm_table/result/rRNA_pm_StaRRA_human_4233_4263.count.xml"
  taskReadFile = "/scratch/stein_lab/shengq2/20200226_4233_4263_michelle_smallRNA_human_v5_byTiger/data_visualization/reads_in_tasks/result/StaRRA_human_4233_4263.NonParallel.TaskReads.csv"
  outputFile="/scratch/stein_lab/shengq2/20200226_4233_4263_michelle_smallRNA_human_v5_byTiger/data_visualization/bacteria_count/result/StaRRA_human_4233_4263.tsv"
else:
  parser = argparse.ArgumentParser(description="Generate smallRNA count from count xml.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-g', '--genomeListFile', action='store', nargs='?', help='Input bacteria genome count xml list file', required=NOT_DEBUG)
  parser.add_argument('-d', '--databaseFile', action='store', nargs='?', help="Original rRNA database count xml file", required=NOT_DEBUG)
  parser.add_argument('-t', '--taskReadFile', action='store', nargs='?', help="Task read count file", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output count file", required=NOT_DEBUG)

  args = parser.parse_args()
  
  print(args)
  
  genomeListFile = args.genomeListFile
  databaseFile = args.databaseFile
  taskReadFile = args.taskReadFile
  outputFile = args.output

logger = logging.getLogger('getBacteriaCount')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def readFileList(fileName):
  result = []
  with open(fileName) as fh:
    for line in fh:
      filepath = line.strip().split('\t', 1)[0]
      result.append(filepath)
  return(result)
   
genomeFiles = readFileList(genomeListFile)

sample_names = set()
result = {}
for genomeFile in genomeFiles:
  logger.info("Parsing " + genomeFile)
  queryMap = {}

  with open(genomeFile, "rt") as fin:
    headers = fin.readline().rstrip().split('\t')
    sample_names.update(headers[1:])
    for line in fin:
      parts = line.rstrip().split('\t')
      query_seq = parts[0]
      for si in range(1, len(parts)):
        sample_name = headers[si]
        query_count = int(parts[si])
        result.setdefault(query_seq, {})[sample_name] = query_count
    
logger.info("Parsing " + databaseFile)
tree = ET.parse(databaseFile)
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
    sample_name = query.get("sample")
    sample_names.add(sample_name)
    result.setdefault(query_seq, {})[sample_name] = query_count

samples = sorted(sample_names)
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
target_r = os.path.basename(rscript)
with open(target_r, "wt") as fout:
  fout.write("outFile='%s'\n" % summaryFile)
  fout.write("parFile1='%s'\n" % summaryFile)
  fout.write("parFile2='%s'\n" % taskReadFile)
  fout.write("setwd('%s')\n\n" % os.path.dirname(os.path.realpath(outputFile)))
  with open(rscript, "rt") as fin:
    for line in fin:
      fout.write(line)

subprocess.call("R --vanilla -f " + target_r, shell=True)

logger.info("done.")
