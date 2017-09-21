import sys
import os
import re
import logging
import argparse
import operator
import xml.etree.ElementTree as ET
from QueryUtils import readQueries

DEBUG = len(sys.argv) == 1

if DEBUG:
  inputFiles ="/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/preprocessing/identical/result/Urine_WT_14_clipped_identical.fastq.dupcount,/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/preprocessing/cutadapt/result/Urine_WT_14_clipped.fastq.short.gz"
  mappedFiles = "/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_genome/bowtie1_bacteria_group1_pm_count/result/Urine_WT_14/Urine_WT_14.bam.count.mapped.xml,/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_genome/bowtie1_bacteria_group2_pm_count/result/Urine_WT_14/Urine_WT_14.bam.count.mapped.xml,/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_genome/bowtie1_fungus_group4_pm_count/result/Urine_WT_14/Urine_WT_14.bam.count.mapped.xml,/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_library/bowtie1_tRNA_pm_count/result/Urine_WT_14/Urine_WT_14.bam.count.mapped.xml,/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_library/bowtie1_rRNA_pm_count/result/Urine_WT_14/Urine_WT_14.bam.count.mapped.xml"
  mappedNames = "bacteria_group1,bacteria_group2,fungus_group4,tRNA,rRNA"
  minCount = 2
  outputFile="/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/unmapped.tsv"
else:
  parser = argparse.ArgumentParser(description="Find the source of unmapped reads.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input identical fastq file')
  parser.add_argument('-b', '--mapped_files', action='store', nargs='?', help="Mapped files, can be either .dupcount or .count.xml, seprated by ','")
  parser.add_argument('-n', '--mapped_names', action='store', nargs='?', help="Mapped names, seprated by ','")
  parser.add_argument('-m', '--min_count', action='store', nargs='?', type=int, default=2, help="Minimum count of read")
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file")

  args = parser.parse_args()
  
  print(args)
  
  inputFiles = args.input
  mappedFiles = args.mapped_files
  mappedNames = args.mapped_names
  minCount = args.min_count
  outputFile= args.output

logger = logging.getLogger('annoateUnmappedReads')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

inputFiles = inputFiles.split(',')
mappedFiles = mappedFiles.split(',')
mappedNames = mappedNames.split(',')

if len(mappedFiles) != len(mappedNames):
  raise Exception("Number of file %d should be same as number of name %d" %(len(mappedFiles), len(mappedNames)))

for f in inputFiles:
  if not os.path.isfile(f):
    raise Exception("File not exists: %s" % f)
  
for f in mappedFiles:
  if not os.path.isfile(f):
    raise Exception("File not exists: %s" % f)

categoryMap = {}
for idx, mappedFile in enumerate(mappedFiles):
  categoryName = mappedNames[idx]
  logger.info("reading %s ..." % mappedFile)
  sequences = readQueries(mappedFile, minCount)
  for query in sequences:
    query.QueryName = re.sub(":CLIP.*", "", query.QueryName)
  categoryMap[categoryName] = sequences
  logger.info("%s : %d" % (categoryName, len(sequences)))

mappedQueries = set()
for sequences in categoryMap.values():
  for item in sequences:
    mappedQueries.add(item.QueryName)


def acceptCandidate(queryItem, mappedQueries):
  seqLen = len(queryItem.Sequence)
  if seqLen < 12 or seqLen > 19:
    return False
  if queryItem.QueryName in mappedQueries:
    return False
  return True

candidateReads = []
for inputFile in inputFiles:
  logger.info("reading %s ..." % inputFile)
  queries = readQueries(inputFile, minCount)
  unmappedQueries = [v for v in queries if acceptCandidate(v, mappedQueries)]
  logger.info("unmapped : %d" % len(unmappedQueries))
  candidateReads.extend(unmappedQueries)

candidateReads.sort(key=operator.attrgetter('QueryCount'), reverse=True)
logger.info("total unmapped : %s" % len(candidateReads))

tempFile = outputFile + ".tmp"
with open(tempFile, "w") as sw:
  sw.write("Sequence\tCount\tLength\t%s\n" % '\t'.join(mappedNames))
  totalCount = len(candidateReads)
  iCount = 0
  for read in candidateReads:
    sequence = read.Sequence
    readCount = read.QueryCount
      
    iCount = iCount + 1
    if iCount % 1000 == 0:
      logger.info("%d / %d" %(iCount, totalCount))
          
    sw.write("%s\t%d\t%d" % (sequence, readCount, len(sequence)))
    for name in mappedNames:
      #logger.info("%s : %s" % (sequence, name))
      sequences = categoryMap[name]
      sw.write("\t")
      for seq in sequences:
        if sequence in seq.Sequence:
          sw.write("%s:%d" % (seq.Sequence, seq.QueryCount))
          break
    sw.write("\n")

if os.path.isfile(outputFile):
  os.remove(outputFile)
os.rename(tempFile, outputFile)

logger.info("done")