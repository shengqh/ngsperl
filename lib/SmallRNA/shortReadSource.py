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
  def __init__(self, sequence, category):
    self.Sequence = sequence
    self.Categories = [category]

def getValue(value):
  return value.TotalCount

def getFilename(value):
  return value[1]

def match(logger, input, names, annotated, output):
  cnames = names.split(",")
  
  logger.info("Reading annotated reads:" + annotated + " ...")
  annotatedReadMap = {}

  with open(annotated, 'r') as annolist:
    iIndex = -1
    for row in annolist:
      parts = row.split('\t')
      annofile = parts[0]
      iIndex = iIndex + 1
      category = cnames[iIndex]
      seqs = []
     
      #icount = 0 
      logger.info("  Reading " + annofile + " ...")
      with open(annofile, 'r') as sr:
        annotatedFiles = sr.readline().rstrip().split('\t')[1:]
        for line in sr:
          parts = line.split('\t', 1)
          seq = parts[0]
          seqs.append(seq)
          #icount = icount + 1
          #if icount >= 1000:
            #break

      annotatedReadMap[category] = seqs

  with open(output, "w") as sw:
    sw.write("Sample\tTotalShortReads\t%s\tUnannotated\n" % ("\t".join(cnames)));

    shortFileList = []
    with open(input, 'r') as sr:
      for line in sr:
        parts = line.rstrip().split('\t')
        shortFileList.append(parts)

    shortFileList = sorted(shortFileList, key=getFilename)
    for parts in shortFileList:
      sampleFile = parts[0]
      sample = parts[1]

      shortReadMap = {}
      for cname in cnames:
        shortReadMap[cname] = 0
      totalCount = 0
      unannoCount = 0

      logger.info("  Reading " + sampleFile + " ...")
      with open(sampleFile, 'r') as fin:
        fin.readline()
        for line in fin:
          reads = line.rstrip().split('\t')
          count = int(reads[1])
          shortSeq = reads[2]
          totalCount = totalCount + count

          bFound = False
          for cname in cnames:
            annotatedReads = annotatedReadMap[cname]

            for annoSeq in annotatedReads:
              if shortSeq in annoSeq:
                bFound = True
                shortReadMap[cname] = shortReadMap[cname] + count
                break

          if not bFound:
            unannoCount = unannoCount + count

      sw.write("%s\t%d\t%s\t%d\n" % (sample, totalCount, "\t".join([str(shortReadMap[cname]) for cname in cnames]), unannoCount))

  logger.info("Done.")

def main():
  parser = argparse.ArgumentParser(description="Matching short reads with annotated reads.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG=False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input short reads', required=NOT_DEBUG)
  parser.add_argument('-a', '--annotated', action='store', nargs='?', help='Input annotated reads', required=NOT_DEBUG)
  parser.add_argument('-n', '--names', action='store', nargs='?', help='Input annotated reads categories, split by ''', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output prefix of matched reads file", required=NOT_DEBUG)
  
  if NOT_DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/cqs/shengq2/vickers/20181205_smallRNA_269_933_2002_human_v4_hg19/data_visualization/short_reads_source_bar/result/human_269_933_2002__fileList1.list"
    args.annotated = "/scratch/cqs/shengq2/vickers/20181205_smallRNA_269_933_2002_human_v4_hg19/data_visualization/short_reads_source_bar/result/human_269_933_2002__fileList2.list"
    args.names = "host smallRNA,host genome,bacteria_group1,bacteria_group2,fungus_group4,virus_group6,nonhost tRNA,nonhost rRNA"
    #args.names = "Host miRNA,Host tRNA"
    args.output = "/scratch/cqs/shengq2/vickers/20181205_smallRNA_269_933_2002_human_v4_hg19/data_visualization/short_reads_source_bar/result/human_269_933_2002.tsv"
  
  logger = logging.getLogger('shortReadSource')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  match(logger, args.input, args.names, args.annotated, args.output)
  
if __name__ == "__main__":
    main()
