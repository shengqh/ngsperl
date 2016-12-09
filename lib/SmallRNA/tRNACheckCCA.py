import sys
import gzip
import io

#inputfile = sys.argv[1]
#originalfile = sys.argv[2]
#outputfile = sys.argv[3]

inputfile="Z:/Shared/Labs/Vickers Lab/Tiger/projects/20150930_TGIRT_tRNA_human/identical/result/KCVH01_clipped_identical.fastq.gz"
originalfile="Z:/Shared/Labs/Vickers Lab/Tiger/data/20150515_tRNA/KCVH1_S6_R1_001.fastq.gz"
outputfile="H:/temp/test_cca.tsv"

ccs={}
if(inputfile.endswith(".gz")):
  f = gzip.open(inputfile, 'rt')
else:
  f = open(inputfile, 'r')

try:
  readCount = 0
  for header in f:
    if not header.startswith("@"):
      continue
    seq = f.__next__().strip()
    ignore = f.__next__().strip()
    score = f.__next__().strip()

    readCount = readCount+1
    if readCount % 10000 == 0:
      print("%d/%d reads end with CC found" % (len(ccs), readCount))

    if seq.endswith("CC"):
      name = header.split(' ')[0]
      ccs[name] = seq
finally:
  f.close()

if(originalfile.endswith(".gz")):
  f = gzip.open(originalfile, 'rt')
else:
  f = open(originalfile, 'r')

try:
  with open(outputfile, "w") as sw:
    readCount = 0
    ccCount = 0
    for header in f:
      if not header.startswith("@"):
        continue
      seq = f.__next__()
      f.__next__()
      f.__next__()

      readCount = readCount + 1
      if readCount % 10000 == 0:
        print("%d/%d reads end with CC processed" % (ccCount, readCount))

      name = header.split(' ')[0]
      sequence = ccs.pop(name, None)

      if sequence == None:
        continue

      ccCount = ccCount + 1

      pos = seq.find(sequence)
      if pos == -1:
        raise ValueError("Cannot find trimmed sequence %s in untrimmed sequence %s of read %s" % (sequence, seq, name))

      if seq[len(sequence)] == 'A':
        sw.write(name + "\n")
finally:
  f.close()

if len(ccs) > 0:
  unfoundFile = outputfile + ".unfound"
  with open(unfoundFile, "w") as fw:
    for key in ccs:
      fw.write(key + "\n")
  raise ValueError("Couldn't find %d reads in untrimmed file, saved to %s" %(len(ccs), unfoundFile))
