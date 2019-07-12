import argparse
import sys
import logging
import os
import csv

class ReadItem:
  def __init__(self, sequence, totalCount):
    self.Sequence = sequence
    self.TotalCount = totalCount
    self.SampleMap = {}

class AnnotationItem:
  def __init__(self, sequence, totalCount, category, counts):
    self.Sequence = sequence
    self.TotalCount = totalCount
    self.Categories = [category]
    self.Counts = counts

def getValue(value):
  return value.TotalCount

def getFilename(value):
  return value[1]

def match(logger, input, names, annotated, maxMapped, maxNumber, minReadCount, minSampleCount, outputPrefix):
  logger.info("Reading short reads:" + input + " ...")
  shortReadMap = {}
  shortReadFiles = []

  shortFileList = []
  with open(input, 'r') as sr:
    for line in sr:
      parts = line.rstrip().split('\t')
      shortFileList.append(parts)
  
  shortFileList = sorted(shortFileList, key=getFilename)
  for parts in shortFileList:
    sampleFile = parts[0]
    sample = parts[1]
    shortReadFiles.append(sample)
    logger.info("  Reading " + sampleFile + " ...")
    with open(sampleFile, 'r') as fin:
      fin.readline()
      for line in fin:
        reads = line.rstrip().split('\t')
        count = int(reads[1])
        seq = reads[2].rstrip()
        if not seq in shortReadMap:
          ri = ReadItem(seq, count)
          shortReadMap[seq] = ri
        else:
          ri = shortReadMap[seq]
          ri.TotalCount += count
        ri.SampleMap[sample] = count
  
  if minSampleCount > 1 or minReadCount > 1:
    shortReads = []
    for read in shortReadMap.values():
      validSampleCount = len([v for v in read.SampleMap.values() if v >= minReadCount])
      if validSampleCount >= minSampleCount:
        shortReads.append(read)
  else:
    shortReads = shortReadMap.values()
  shortReads = sorted(shortReads, key=getValue, reverse=True)
  
  if len(shortReads) > maxNumber:
    shortReads = shortReads[0:maxNumber]

  logger.info("Reading max mapped reads:" + maxMapped + " ...")
  maxmappedReads = {}
  with open(maxMapped, 'r') as sr:
    for line in sr:
      parts = line.split('\t')
      logger.info("  Reading " + parts[0] + " ...")
      with open(parts[0], 'r') as fin:
        while True:
          qname = fin.readline().rstrip()
          if not qname:
            break
          seq = fin.readline()
          fin.readline()
          fin.readline()
  
          if qname.endswith("_"):
            maxmappedReads[seq.rstrip()] = 1
        
  cnames = names.split(",")
  
  logger.info("Reading annotated reads:" + annotated + " ...")
  annotatedReadMap = {}
  annotatedFiles = []
  with open(annotated, 'r') as annolist:
    iIndex = -1
    for row in annolist:
      parts = row.split('\t')
      annofile = parts[0]
      iIndex = iIndex + 1
      category = cnames[iIndex]
      
      logger.info("  Reading " + annofile + " ...")
      with open(annofile, 'r') as sr:
        annotatedFiles = sr.readline().rstrip().split('\t')[1:]
        for line in sr:
          parts = line.rstrip().split('\t')
          seq = parts[0]
          if seq not in annotatedReadMap:
            totalCount = sum(int(p) for p in parts[1:])
            annotatedReadMap[seq] = AnnotationItem(seq, totalCount, category, parts[1:])
          else:
            annotatedReadMap[seq].Categories.append(category)
            
  annotatedReads = sorted(annotatedReadMap.values(), key=getValue, reverse=True)
  
  output = outputPrefix + ".tsv"
  logger.info("Writing explain result:" + output + " ...")
  with open(output, "w") as sw:
    sw.write("ShortRead\tShortReadCount\tShortReadLength\t" + "\t".join(["SRS_" + f for f in shortReadFiles]) + "\tIsMaxMapped\tParentRead\tParentReadCount\tParentReadCategory\t" + "\t".join(["PRS_" + f for f in annotatedFiles]) + "\n")
    emptyAnnotation = "\t\t\t\t" + "\t".join(["" for af in annotatedFiles]) + "\n"
    for shortRead in shortReads:
      shortSeq = shortRead.Sequence
      shortSeqCount = shortRead.TotalCount
      seqMap = shortRead.SampleMap
      sw.write("%s\t%s\t%d" % (shortSeq, shortSeqCount, len(shortSeq)))
      for fname in shortReadFiles:
        if fname in seqMap:
          sw.write("\t%s" % seqMap[fname])
        else:
          sw.write("\t0")
      
      sw.write("\t" + str(shortSeq in maxmappedReads))
          
      bFound = False
      for annotatedRead in annotatedReads:
        annoSeq = annotatedRead.Sequence
        if shortSeq in annoSeq:
          bFound = True
          sw.write("\t%s\t%s\t%s\t%s\n" % (annoSeq, annotatedRead.TotalCount, "/".join(annotatedRead.Categories[0]), "\t".join(annotatedRead.Counts)))
          break
      
      if not bFound:
        sw.write(emptyAnnotation)
          
  logger.info("Done.")

def main():
  parser = argparse.ArgumentParser(description="Matching short reads with annotated reads.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG=False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input short reads', required=NOT_DEBUG)
  parser.add_argument('-m', '--maxMapped', action='store', nargs='?', help='Input reads exceed maximum mapping to genome', required=NOT_DEBUG)
  parser.add_argument('-a', '--annotated', action='store', nargs='?', help='Input annotated reads', required=NOT_DEBUG)
  parser.add_argument('-n', '--names', action='store', nargs='?', help='Input annotated reads categories, split by ''', required=NOT_DEBUG)
  parser.add_argument('--maxNumber', action='store', default=100, nargs='?', help='Input number of top short reads for annotation')
  parser.add_argument('--minReadCount', action='store', default=3, nargs='?', help='Input minimum copy of short reads in sample for annotation')
  parser.add_argument('--minSampleCount', action='store', default=2, nargs='?', help='Input minimum number of sample with valid read count')
  parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output prefix of matched reads file", required=NOT_DEBUG)
  
  if NOT_DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/data_visualization/short_reads_source/result/match__fileList1.list"
    args.maxMapped = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/data_visualization/short_reads_source/result/match__fileList2.list"
    args.annotated = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/data_visualization/short_reads_source/result/match__fileList3.list"
    args.names = "Host miRNA,Host tRNA,Host snRNA,Host snoRNA,Host rRNA,Host other small RNA,Host Genome,Microbiome Bacteria,Environment Bacteria,Fungus,Non host tRNA,Non host rRNA"
    #args.names = "Host miRNA,Host tRNA"
    args.output = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/data_visualization/short_reads_source/result/match2"
  
  logger = logging.getLogger('explainShortReads')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  match(logger, args.input, args.names, args.annotated, args.maxMapped, args.maxNumber, args.minReadCount, args.minSampleCount, args.output)
  
if __name__ == "__main__":
    main()
