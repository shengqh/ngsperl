import sys
import os
import logging
import argparse
import pysam
import subprocess
from DupCountUtils import readDupCountMap
from FastqUtils import readFastqQueryNames
from FileListUtils import readFileMap

DEBUG = 0

if DEBUG:
  inputFile="/scratch/cqs/shengq2/vickers/20180410_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_genome/bowtie1_bacteria_group1_pm_mappedreads_host_mismatch_table/result/KCV_3018_77_78_79__fileList1.list"
  fastqFile="/scratch/cqs/shengq2/vickers/20180410_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_genome/bowtie1_bacteria_group1_pm_mappedreads_host_mismatch_table/result/KCV_3018_77_78_79__fileList2.list"
  countFile="/scratch/cqs/shengq2/vickers/20180410_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_genome/bowtie1_bacteria_group1_pm_mappedreads_host_mismatch_table/result/KCV_3018_77_78_79__fileList3.list"
  outputFile="/scratch/cqs/shengq2/temp/hostgenome.tsv"
  maxMismatch=2
else:
  parser = argparse.ArgumentParser(description="Extract mapped reads to Fastq.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bam file list')
  parser.add_argument('-f', '--fastq', action='store', nargs='?', help="Original fastq file list ")
  parser.add_argument('-c', '--count', action='store', nargs='?', help="Original dupcount file list")
  parser.add_argument('-m', '--maxMismatch', action='store', nargs='?', default=2, type=int, help="Maximum number of mismatch")
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output result file")

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  fastqFile = args.fastq
  countFile = args.count
  outputFile = args.output
  maxMismatch = args.maxMismatch

logger = logging.getLogger('bamMismatchCount')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

bamFiles = readFileMap(inputFile)
fastqFiles = readFileMap(fastqFile)
countFiles = readFileMap(countFile)

mismatchRange = range(0, maxMismatch+1)
misMap = {}
totalMap = {}
samples = bamFiles.keys()
for sample in samples:
  bamFile = bamFiles[sample]
  fastqFile = fastqFiles[sample]
  countFile = countFiles[sample]
  counted = set()
  
  logger.info("Reading counts from " + countFile)
  counts = readDupCountMap(countFile)

  logger.info("Reading original reads from " + fastqFile)
  queryNames = readFastqQueryNames(fastqFile)

  totalCount = sum(counts[q] for q in queryNames)
  totalMap[sample] = totalCount

  misCount = {}
  misMap[sample] = misCount
  for mis in mismatchRange:
    misCount[mis] = 0

  logger.info("Iterating BAM file " + bamFile)
  with pysam.Samfile(bamFile) as sam:
    for read in sam.fetch():
      if read.query_name in counted:
        continue
      counted.add(read.query_name)
      
      numberOfMismatch = 0
      for tag in read.tags:
        if tag[0] == "NM":
          numberOfMismatch = tag[1]
          break
      
      qcount = counts[read.query_name]
      if numberOfMismatch in misCount:
        newCount = misCount[numberOfMismatch] + qcount
        misCount[numberOfMismatch] = newCount
      else:
        misCount[numberOfMismatch] = qcount
  
with open(outputFile, "w") as sw:
  sw.write("Category\t%s\n" % "\t".join(samples))
  for mis in mismatchRange:
    sw.write("Mismatch_%d\t%s\n" % (mis, "\t".join([str(misMap[sample][mis]) for sample in samples])))
  sw.write("Unmapped\t%s\n" % "\t".join([str(totalMap[sample] - sum(misMap[sample].values())) for sample in samples]))

rscript = os.path.realpath(__file__) + ".R"
subprocess.call ("R --vanilla -f " + rscript + " --args \"" + outputFile + "\" \"" + outputFile + ".pdf\"", shell=True)

logger.info("Done")
