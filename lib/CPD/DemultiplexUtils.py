import gzip

def readFileMap(fileName):
  result = {}
  with open(fileName) as fh:
    for line in fh:
      filepath, name = line.strip().split('\t', 1)
      result[name] = filepath.strip()
  return(result)
 
def demultiplex(inputFile, outputFilePrefix, configFile, args, logger):
  sampleSeqMap = readFileMap(configFile)
  seqSampleMap = {sampleSeqMap[k]:k for k in sampleSeqMap.keys()}
  
  print(seqSampleMap)
  seqFileMap = {}
  barcodeLength = 0
  for seq in seqSampleMap.keys():
    barcodeLength = len(seq)
    seqFile = outputFilePrefix + "." + seqSampleMap[seq] + ".fastq.gz"
    seqFileMap[seq] = gzip.open(seqFile, 'wt')
  
  logger.info("reading input file: " + inputFile + " ...")
  count = 0
  with open(inputFile, "r") as fin:
    while(True):
      query = fin.readline()
      if not query:
        break
      
      seq = fin.readline().rstrip()
      skipline = fin.readline()
      score = fin.readline().rstrip()
      
      count = count + 1
      if count % 100000 == 0:
        logger.info("%d processed" % count)

      barcode = seq[0:barcodeLength]
      if barcode in seqFileMap:
        newseq = seq[barcodeLength:-1]
        newscore = score[barcodeLength:-1]
        fout = seqFileMap[barcode]
        fout.write(query)
        fout.write(newseq + '\n')
        fout.write(skipline)
        fout.write(newscore + '\n')
        
  for seq in seqFileMap.keys():
    seqFileMap[seq].close()
  
  logger.info("demultiplex done.")
