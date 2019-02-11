import argparse
import sys
import logging
import os
import csv

def getKey(value):
  return int(value[1])

def match(input, annotated, output, logger):
  shortReads = []
  shortHeader = ""
  with open(input, 'r') as sr:
    shortHeader = sr.readline().rstrip()
    for line in sr:
      parts = line.rstrip().split('\t')
      if int(parts[1]) > 1:
        shortReads.append(parts)
  
  shortReads.sort(key=getKey, reverse=True)
  
  annotatedReads = []
  annotatedHeader = ""
  with open(annotated, 'r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    columns = csv_reader.next()
    columns[0] = "AnnotatedRead"
    annotatedHeader = '\t'.join(columns)
    #  print(annotatedHeader)
    for row in csv_reader:
      annotatedReads.append(row)
  
  with open(output, "w") as sw:
    sw.write(shortHeader + "\t" + annotatedHeader + "\n")
    for shortRead in shortReads:
      shortSeq = shortRead[2]
      for annotatedRead in annotatedReads:
        annoSeq = annotatedRead[0]
        if shortSeq in annoSeq:
          sw.write('\t'.join(shortRead) + '\t' + '\t'.join(annotatedRead) + '\n')        
  
  

def main():
  parser = argparse.ArgumentParser(description="Matching short reads with annotated reads.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG=False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input short reads', required=NOT_DEBUG)
  parser.add_argument('-a', '--annotated', action='store', nargs='?', help='Input annotated reads', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output matched reads file", required=NOT_DEBUG)
  
  if NOT_DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/host_genome/bowtie1_genome_unmapped_reads/result/cell_5mM_1_clipped_identical.short.fastq.dupcount";
    args.annotated = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/data_visualization/sequence_mapped_in_categories/result/human_269_933_2002.ReadsMapping.Summary.csv";
    args.output = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/host_genome/bowtie1_genome_unmapped_reads/result/cell_5mM_1_clipped_identical.short.fastq.dupcount.matched.tsv"
  
  logger = logging.getLogger('explainShortReads')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  match(args.input, args.annotated, args.output, logger)
  
if __name__ == "__main__":
    main()
