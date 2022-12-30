import sys
import gzip
import os

DEBUG = 0

if DEBUG:
  inputFile="Z:/Shared/Labs/Vickers Lab/Tiger/projects/20150930_TGIRT_tRNA_human/identical/result/KCVH01_clipped_identical.fastq.gz"
  originalFile="Z:/Shared/Labs/Vickers Lab/Tiger/data/20150515_tRNA/KCVH1_S6_R1_001.fastq.gz"
  outputFile="H:/temp/test_cca.tsv"
else:
  inputFile = sys.argv[1]
  originalFile = sys.argv[2]
  outputFile = sys.argv[3]

ccs={}
if(inputFile.endswith(".gz")):
  f = gzip.open(inputFile, 'rt')
else:
  f = open(inputFile, 'r')

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

    readCount = readCount+1
    if readCount % 10000 == 0:
      print("%d/%d reads end with CC found" % (len(ccs), readCount))

    if seq.endswith("CC"):
      name = header.split(' ')[0]
      ccs[name] = seq
finally:
  f.close()

if(originalFile.endswith(".gz")):
  f = gzip.open(originalFile, 'rt')
else:
  f = open(originalFile, 'r')

try:
  tempFile = outputFile + ".tmp"
  with open(tempFile, "w") as sw:
    readCount = 0
    ccCount = 0
    while True:
      header = f.readline()
      if '' == header:
        break

      if not header.startswith("@"):
        continue
      seq = f.readline()
      f.readline()
      f.readline()

      readCount = readCount + 1
      if readCount % 100000 == 0:
        print("%d/%d reads end with CC processed" % (ccCount, readCount))

      name = header.split(' ')[0]
      sequence = ccs.pop(name, None)

      if sequence is None:
        continue

      ccCount = ccCount + 1

      pos = seq.find(sequence)
      if pos == -1:
        raise ValueError("Cannot find trimmed sequence %s in untrimmed sequence %s of read %s" % (sequence, seq, name))

      if seq[len(sequence)] == 'A':
        sw.write(name + "\n")

  if os.path.isfile(outputFile):
    os.remove(outputFile)
  os.rename(tempFile, outputFile)
finally:
  f.close()

if len(ccs) > 0:
  unfoundFile = outputFile + ".unfound"
  with open(unfoundFile, "w") as fw:
    for key in ccs:
      fw.write(key + "\n")
  raise ValueError("Couldn't find %d reads in untrimmed file, saved to %s" %(len(ccs), unfoundFile))
