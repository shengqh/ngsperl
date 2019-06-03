import argparse
import sys
import logging
import os
import csv

def getTotalCounts(value):
  return value[len(value)-1]

def readCountTable(fileName):
  result = []
  headers = []
  with open(fileName, 'r') as sr:
    headers = sr.readline().rstrip().split('\t')
    for line in sr:
      parts = line.rstrip().split('\t')
      result.append(parts)
  return [headers, result]    
  
def match(logger, input, namesStr, annotated, maxmapped, output):
  logger.info("Reading max mapped reads:" + maxmapped + " ...")
  
  names = namesStr.split(',')
  maxmappedReads = {}
  with open(maxmapped, 'r') as srfile:
    for fileline in srfile:
      parts = fileline.split('\t')
      with open(parts[0], 'r') as sr:
        lines = []
        for line in sr:
          lines.append(line)
          if len(lines) == 4:
            maxmappedReads[lines[1].rstrip()] = 1
            lines = []
  
  logger.info("Reading short reads:" + input + " ...")
  shortReadPair = readCountTable(input)
  shortHeaders = shortReadPair[0]
  shortReads = shortReadPair[1]
  
  #keep the read found in at least two samples
  shortReads = [ parts for parts in shortReads if sum(parts[index] != "0" for index in range(1, len(parts))) > 1 ]
  
  shortHeaders.append("IsMaxMapped")
  shortHeaders.append("TotalReads")
  for parts in shortReads:
    counts = sum(int(parts[index]) for index in range(1, len(parts)))
    parts.append(str(parts[0] in maxmappedReads))
    parts.append(counts)
    
  shortReads.sort(key=getTotalCounts, reverse=True)
  
  logger.info("Reading annotated reads:" + annotated + " ...")
  annotatedReads = {}
  annotatedHeader = "\t".join(names)
  with open(annotated, 'r') as annoFileList:
    for name in names:
      annoFile = annoFileList.readline().split('\t')[0]
      annotatedReads[name] = readCountTable(annoFile)
  
  logger.info("Writing explain result:" + output + " ...")
  with open(output, "w") as sw:
    sw.write(shortHeader + "\t" + annotatedHeader + "\n")
    for shortRead in shortReads:
      shortSeq = shortRead[0]
      bFound = False
      for annotatedRead in annotatedReads:
        annoSeq = annotatedRead[0]
        if shortSeq in annoSeq:
          bFound = True
          sw.write('\t'.join(shortRead) + '\t' + '\t'.join(annotatedRead) + '\n')
      
      if not bFound:
          sw.write('\t'.join(shortRead) + '\t' + emptyAnnotation + '\n')
          
  logger.info("Done.")

def main():
  parser = argparse.ArgumentParser(description="Matching short reads with annotated reads.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG=True
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input short reads', required=NOT_DEBUG)
  parser.add_argument('-a', '--annotated', action='store', nargs='?', help='Input annotated reads', required=NOT_DEBUG)
  parser.add_argument('-m', '--maxmapped', action='store', nargs='?', help='Input reads exceed maximum mapping to genome', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output matched reads file", required=NOT_DEBUG)
  
  if NOT_DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/cqs/shengq2/vickers/20190304_smallRNA_3065_DM_1_mouse_V4T_tRH/host_genome/bowtie1_genome_host_too_short_reads_table/result/DM_3065_mouse.count"
    args.annotated = "/scratch/cqs/shengq2/vickers/20190304_smallRNA_3065_DM_1_mouse_V4T_tRH/data_visualization/short_reads_source/result/DM_3065_mouse__fileList3.list"
    args.maxmapped = "/scratch/cqs/shengq2/vickers/20190304_smallRNA_3065_DM_1_mouse_V4T_tRH/data_visualization/short_reads_source/result/DM_3065_mouse__fileList2.list"
    args.output = "/scratch/cqs/shengq2/vickers/20190304_smallRNA_3065_DM_1_mouse_V4T_tRH/data_visualization/short_reads_source/result/DM_3065_mouse.tsv";
    args.names="Host miRNA,Host tRNA,Host snRNA,Host snoRNA,Host rRNA,Host other small RNA,Host Genome,Microbiome Bacteria,Environment Bacteria,Fungus,Non host tRNA,Non host rRNA"
#     args.annotated = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/data_visualization/sequence_mapped_in_categories/result/human_269_933_2002.ReadsMapping.Summary.csv";
#     args.maxmapped = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/host_genome/bowtie1_genome_1mm_NTA/result/cell_5mM_1.bam.max.txt";
#     args.output = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/host_genome/bowtie1_genome_unmapped_reads/result/cell_5mM_1_clipped_identical.short.fastq.dupcount.matched.tsv"
  
  logger = logging.getLogger('explainShortReads')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  match(logger, args.input, args.names, args.annotated, args.maxmapped, args.output)
  
if __name__ == "__main__":
    main()
