import sys
import gzip
import os
import logging

logger = logging.getLogger('fastqSmallRnaNTA')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

DEBUG = 1
NTA_TAG = ":CLIP_"

if DEBUG:
  inputFile="Z:/Shared/Labs/Vickers Lab/Tiger/projects/20150930_TGIRT_tRNA_human/identical/result/KCVH01_clipped_identical.fastq.gz"
  outputFile="H:/temp/KCVH01_clipped_identical_NTA.fastq.gz"
  minReadLength=16
else:
  inputFile = sys.argv[1]
  outputFile = sys.argv[2]
  minReadLength = int(sys.argv[3])

def readCountMap(countFile, logger):
  result = {}
  if countFile != None:
    logger.info("Reading count from %s ..." % countFile)
    readCount = 0
    with open(countFile, "r") as f:
      header = f.readline()
      for line in f:
        readCount = readCount + 1
        if readCount % 10000 == 0:
          logger.info("%d read count readed" % readCount)
        parts = line.rstrip().split('\t')
        result["@" + parts[0]] = parts[1]
  return result

logger.info("Processing queries in %s ..." % inputFile)
if inputFile.endswith(".gz"):
  f = gzip.open(inputFile, 'rt')
else:
  f = open(inputFile, 'r')

try:
  tmpFile=outputFile + ".tmp"
  if outputFile.endswith(".gz"):
    fw = gzip.open(tmpFile, "wt")
  else:
    fw = open(tmpFile, "w")

  try:
    readCount = 0
    rng3 = range(0,4)
    rng4 = range(0,5)
    while True:
      header = f.readline()
      if '' == header:
        break

      if not header.startswith("@"):
        continue

      seq = f.readline().strip()
      ignore = f.readline().strip()
      score = f.readline().strip()
      seqLength = len(seq)

      readCount = readCount + 1
      if readCount % 10000 == 0:
        logger.info("%d reads processed" % readCount)

      qname = header.split(' ')[0]

      curRange = rng4 if seq.endswith("CCAA") else rng3
      for i in curRange:
        newLength = seqLength - i
        if newLength < minReadLength:
          break

        clipped = "" if i == 0 else seq[newLength]
        header = qname + NTA_TAG + clipped

        fw.write(header + "\n")
        fw.write(ignore + "\n")
        fw.write(seq[0:newLength] + "\n")
        fw.write(score[0:newLength] + "\n")
    fw.close()
    if os.path.isfile(outputFile):
      os.remove(outputFile)
      os.rename(tmpFile, outputFile)
    logger.info("Generating NTA reads successed!")
  except BaseException as e:
    logger.error('Failed to generate NTA reads: %s' % str(e))
    fw.close()
finally:
  f.close()
