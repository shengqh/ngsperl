import sys
import gzip
import os
import logging
import argparse

DEBUG = 0

if DEBUG:
  input="Z:/Shared/Labs/Vickers Lab/Tiger/projects/20150930_TGIRT_tRNA_human/identical/result/KCVH01_clipped_identical.fastq.gz"
  output="H:/temp/KCVH01_clipped_identical_NTA.fastq.gz"
else:
  parser = argparse.ArgumentParser(description="Generate smallRNA NTA read for Fastq file.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Fastq file')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output Fastq file")

  args = parser.parse_args()
  
  print(args)
  
  input = args.input
  output = args.output

inputFiles = input.split(",")
outputFiles = output.split(",")

logger = logging.getLogger('removeFastqTerminalN')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

for idx in range(0, len(inputFiles)):
  inputFile = inputFiles[idx]
  outputFile = outputFiles[idx]

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

        while(seq.startswith('N')):
          seq = seq[1:]
          score = score[1:]

        while(seq.endswith('N')):
          seq = seq[:-1]
          score = score[:-1]

        fw.write(header.rstrip() + "\n")
        fw.write(seq + "\n")
        fw.write(ignore + "\n")
        fw.write(score + "\n")
      fw.close()
      if os.path.isfile(outputFile):
        os.remove(outputFile)
      os.rename(tmpFile, outputFile)
      logger.info("Remove terminal N successed!")
    except BaseException as e:
      logger.error('Failed to remove terminal N: %s' % str(e))
      fw.close()
  finally:
    f.close()
