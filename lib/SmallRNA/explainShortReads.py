import argparse
import sys
import logging
import os
import csv

def getKey(value):
  return int(value[1])

def match(logger, input, annotated, maxmapped, output):
  logger.info("Reading max mapped reads:" + maxmapped + " ...")
  maxmappedReads = {}
  with open(maxmapped, 'r') as sr:
    lines = []
    for line in sr:
      lines.append(line)
      if len(lines) == 4:
        maxmappedReads[lines[1].rstrip()] = 1
        lines = []
  
  logger.info("Reading short reads:" + maxmapped + " ...")
  shortReads = []
  shortHeader = ""
  with open(input, 'r') as sr:
    shortHeader = sr.readline().rstrip() + "\tIsMaxMapped"
    for line in sr:
      parts = line.rstrip().split('\t')
      if int(parts[1]) > 1:
        parts.append(str(parts[2] in maxmappedReads))
        shortReads.append(parts)
  
  shortReads.sort(key=getKey, reverse=True)
  
  logger.info("Reading annotated reads:" + annotated + " ...")
  annotatedReads = []
  annotatedHeader = ""
  emptyAnnotation = ""
  with open(annotated, 'r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    columns = csv_reader.next()
    columns[0] = "AnnotatedRead"
    annotatedHeader = '\t'.join(columns)
    emptyAnnotation = '\t'.join('' for c in columns) 
    #  print(annotatedHeader)
    for row in csv_reader:
      annotatedReads.append(row)
  
  logger.info("Writing explain result:" + output + " ...")
  with open(output, "w") as sw:
    sw.write(shortHeader + "\t" + annotatedHeader + "\n")
    for shortRead in shortReads:
      shortSeq = shortRead[2]
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
  
  DEBUG=False
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
    args.input = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/host_genome/bowtie1_genome_unmapped_reads/result/cell_5mM_1_clipped_identical.short.fastq.dupcount";
    args.annotated = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/data_visualization/sequence_mapped_in_categories/result/human_269_933_2002.ReadsMapping.Summary.csv";
    args.maxmapped = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/host_genome/bowtie1_genome_1mm_NTA/result/cell_5mM_1.bam.max.txt";
    args.output = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/host_genome/bowtie1_genome_unmapped_reads/result/cell_5mM_1_clipped_identical.short.fastq.dupcount.matched.tsv"
  
  logger = logging.getLogger('explainShortReads')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  match(logger, args.input, args.annotated, args.maxmapped, args.output)
  
if __name__ == "__main__":
    main()
