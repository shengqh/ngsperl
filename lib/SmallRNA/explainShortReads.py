import argparse
import sys
import logging
import os
import csv

def getValue(value):
  return int(value[1][0])

def match(logger, input, names, annotated, maxMapped, maxNumber, output):
  logger.info("Reading short reads:" + input + " ...")
  shortReadMap = {}
  shortReadFiles = []
  with open(input, 'r') as sr:
    for line in sr:
      parts = line.split('\t')
      shortReadFiles.append(parts[1].rstrip())
      logger.info("  Reading " + parts[0] + " ...")
      with open(parts[0], 'r') as fin:
        fin.readline()
        for line in fin:
          reads = line.split('\t')
          seq = reads[2].rstrip()
          if not seq in shortReadMap:
            shortReadMap[seq] = [int(reads[1]), {parts[1]:reads[1]}]
          else:
            curArray = shortReadMap[seq]
            curArray[0] = curArray[0] + int(reads[1])
            curArray[1][parts[1]] = reads[1]
  
  shortReads = sorted(shortReadMap.items(), key=getValue, reverse=True)
  shortReads = shortReads[0:maxNumber]

  logger.info("Reading max mapped reads:" + maxMapped + " ...")
  maxmappedReads = {}
  with open(maxMapped, 'r') as sr:
    for line in sr:
      parts = line.split('\t')
      logger.info("  Reading " + parts[0] + " ...")
      with open(parts[0], 'r') as fin:
        while True:
          qname = fin.readline()
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
      
      logger.info("  Reading " + annofile + " ...")
      with open(annofile, 'r') as sr:
        annotatedFiles = sr.readline().rstrip().split('\t')[1:]
        for line in sr:
          parts = line.split('\t')
          if parts[0] not in annotatedReadMap:
            annotatedReadMap[parts[0]] = [sum(int(p) for p in parts[1:]), [cnames[iIndex]],  parts[1:]]
          else:
            annotatedReadMap[parts[0]][1].append(cnames[iIndex])
            
  annotatedReads = sorted(annotatedReadMap.items(), key=getValue, reverse=True)
  
  logger.info("Writing explain result:" + output + " ...")
  with open(output, "w") as sw:
    sw.write("ShortReads\tShortReadLength\t" + "\t".join(shortReadFiles) + "\tIsMaxMapped\t" + "\t".join(cnames) + "\t" + "\t".join(annotatedFiles) + "\n")
    sw.write(shortHeader + "\t" + annotatedHeader + "\n")
    for shortRead in shortReads:
      shortSeq = shortRead[0]
      seqMap = shortRead[1][1]
      sw.write("%s\t%d" % (shortSeq, len(shortSeq)))
      for fname in shortReadFiles:
        if fname in seqMap:
          sw.write("\t%s" % seqMap[fname])
        else:
          sw.write("\t0")
      
      sw.write("\t" + str(shortSeq in maxmappedReads))
          
      bFound = False
      for annotatedRead in annotatedReads:
        annoSeq = annotatedRead[0]
        if shortSeq in annoSeq:
          bFound = True
          annoInfo = annotatedRead[1]
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
  parser.add_argument('-m', '--maxMapped', action='store', nargs='?', help='Input reads exceed maximum mapping to genome', required=NOT_DEBUG)
  parser.add_argument('-a', '--annotated', action='store', nargs='?', help='Input annotated reads', required=NOT_DEBUG)
  parser.add_argument('-n', '--names', action='store', nargs='?', help='Input annotated reads categories, split by ''', required=NOT_DEBUG)
  parser.add_argument('--maxNumber', action='store', default=500, nargs='?', help='Input number of top short reads for annotation')
  parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output matched reads file", required=NOT_DEBUG)
  
  if NOT_DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/data_visualization/short_reads_source/result/match__fileList1.list";
    args.maxMapped = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/data_visualization/short_reads_source/result/match__fileList2.list";
    args.annotated = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/data_visualization/short_reads_source/result/match__fileList3.list";
    args.names = "Host miRNA,Host tRNA,Host snRNA,Host snoRNA,Host rRNA,Host other small RNA,Host Genome,Microbiome Bacteria,Environment Bacteria,Fungus,Non host tRNA,Non host rRNA"
    args.output = "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180809_smallRNA_269_933_2002_human/data_visualization/short_reads_source/result/match.tsv"
  
  logger = logging.getLogger('explainShortReads')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  match(logger, args.input, args.names, args.annotated, args.maxMapped, args.maxNumber, args.output)
  
if __name__ == "__main__":
    main()
