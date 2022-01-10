import argparse
import sys
import logging
import os
import csv

class ReadItem:
  def __init__(self, sequence):
    self.Sequence = sequence
    self.TotalCount = 0
    self.SampleMap = {}

class AnnotationItem:
  def __init__(self, sequence, category):
    self.Sequence = sequence
    self.TotalCount = 0
    self.Categories = [category]
    self.SampleMap = {}

def getValue(value):
  return value.TotalCount

def getFilename(value):
  return value[1]

def match(logger, input, names, annotated, maxMapped, maxNumber, minReadCount, minSampleCount, outputFile):
  logger.info("Reading short reads:" + input + " ...")
  shortReadMap = {}
  with open(input, 'r') as sr:
    samples = sr.readline().rstrip().split('\t')
    for line in sr:
      parts = line.rstrip().split('\t')
      ri = ReadItem(parts[0])
      shortReadMap[parts[0]] = ri
      for si in range(1, len(samples)):
        sample = samples[si]
        count = int(parts[si])
        ri.SampleMap[sample] = count
        ri.TotalCount += count
  
  samples = samples[1:]

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
  annotatedSamples = []
  with open(annotated, 'r') as annolist:
    iIndex = -1
    for row in annolist:
      parts = row.split('\t')
      annofile = parts[0]
      iIndex = iIndex + 1
      category = cnames[iIndex]
      start_index = 3 if category.startswith('Host ') else 1
      
      logger.info("  Reading " + category + " : " + annofile + " ...")
      with open(annofile, 'r') as sr:
        annotatedSamples = sr.readline().rstrip().split('\t')
        for line in sr:
          parts = line.rstrip().split('\t')
          seq = parts[0]
          if seq not in annotatedReadMap:
            ai = AnnotationItem(seq, category)
            annotatedReadMap[seq] = ai
            for idx in range(start_index, len(annotatedSamples)):
              sample = annotatedSamples[idx]
              count = int(parts[idx])
              ai.SampleMap[sample] = count
              ai.TotalCount += count
          else:
            annotatedReadMap[seq].Categories.append(category)
  annotatedSamples = annotatedSamples[1:]            
  annotatedReads = sorted(annotatedReadMap.values(), key=getValue, reverse=True)
  
  logger.info("Writing explain result:" + outputFile + " ...")
  with open(outputFile, "w") as sw:
    sw.write("ShortRead\tShortReadCount\tShortReadLength\t" + "\t".join(["SRS_" + f for f in samples]) + "\tIsMaxMapped\tParentRead\tParentReadCount\tParentReadCategory\t" + "\t".join(["PRS_" + f for f in annotatedSamples]) + "\n")
    emptyAnnotation = "\t\t\t\t" + "\t".join(["" for af in annotatedSamples]) + "\n"
    for shortRead in shortReads:
      shortSeq = shortRead.Sequence
      shortSeqCount = shortRead.TotalCount
      seqMap = shortRead.SampleMap
      sw.write("%s\t%s\t%d" % (shortSeq, shortSeqCount, len(shortSeq)))
      for fname in samples:
        if fname in seqMap:
          sw.write("\t%s" % seqMap[fname])
        else:
          sw.write("\t0")
      
      sw.write("\t" + str(shortSeq in maxmappedReads))
          
      bFound = False
      for ai in annotatedReads:
        annoSeq = ai.Sequence
        if shortSeq in annoSeq:
          bFound = True
          sw.write("\t%s\t%s\t%s\t%s\n" % (annoSeq, ai.TotalCount, "/".join(ai.Categories), "\t".join([str(ai.SampleMap[sample]) if sample in ai.SampleMap else "0" for sample in annotatedSamples])))
          break
      
      if not bFound:
        sw.write(emptyAnnotation)
          
  logger.info("Done.")

def main():
  parser = argparse.ArgumentParser(description="Matching short reads with annotated reads.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG=False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input short reads table', required=NOT_DEBUG)
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
    args.input = "/scratch/vickers_lab/projects/20210503_6280_RA_smRNAseq_mouse_byTiger/host_genome/bowtie1_genome_short_reads_table/result/RA_6280_mouse.count.txt"
    args.maxMapped = "/scratch/vickers_lab/projects/20210503_6280_RA_smRNAseq_mouse_byTiger/data_visualization/short_reads_source/result/RA_6280_mouse__fileList1.list"
    args.annotated = "/scratch/vickers_lab/projects/20210503_6280_RA_smRNAseq_mouse_byTiger/data_visualization/short_reads_source/result/RA_6280_mouse__fileList2.list"
    args.names = "Host miRNA,Host tRNA,Host snRNA,Host snoRNA,Host rRNA,Host other small RNA,Host Genome,Microbiome Bacteria,Environment Bacteria,Fungus,Algae,Virus,Non host tRNA,Non host rRNA"
    args.output = "/scratch/vickers_lab/projects/20210503_6280_RA_smRNAseq_mouse_byTiger/data_visualization/short_reads_source/result/RA_6280_mouse.tsv"
  
  logger = logging.getLogger('shortReadSource')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  match(logger, args.input, args.names, args.annotated, args.maxMapped, args.maxNumber, args.minReadCount, args.minSampleCount, args.output)
  
if __name__ == "__main__":
    main()
