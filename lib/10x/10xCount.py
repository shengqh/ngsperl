import pysam
import argparse
import sys
import logging
import os
  
def count(outputFile, inputFile, logger):
  gapSampleCountMap = {}
  otherSampleCountMap = {}
  with pysam.AlignmentFile(inputFile, "rb") as samfile:
    processed = 0
    for read in samfile.fetch(until_eof=True):
      processed += 1
      if processed % 1000000 == 0:
        logger.info("processed %d" % (processed))
        
      if read.is_unmapped:
        continue
      
      qname = read.query_name
      
      gap = 0
      for ct in read.cigartuples:
        if ct[0] == 3:
          gap = ct[1]
          break
        
      if gap > 150:
        continue
      
      sample = qname.split(':')[1]
      if gap == 150:
        if sample in gapSampleCountMap:
          gapSampleCountMap[sample] = gapSampleCountMap[sample] + 1
        else:
          gapSampleCountMap[sample] = 1
      else:
        if sample in otherSampleCountMap:
          otherSampleCountMap[sample] = otherSampleCountMap[sample] + 1
        else:
          otherSampleCountMap[sample] = 1
        
  samples = sorted(set(gapSampleCountMap.keys()).union(set(otherSampleCountMap.keys())))
  print(samples)
  with open(outputFile, "wt") as sw:
    sw.write("Sample\tGap\tOther\n")
    for sample in samples:
      sw.write("%s\t%d\t%d\n" % (sample, 
                                 gapSampleCountMap[sample] if sample in gapSampleCountMap else 0,
                                 otherSampleCountMap[sample] if sample in otherSampleCountMap else 0
                                 ))
  
  logger.info("Done.")
        
def main():
  parser = argparse.ArgumentParser(description="Count the reads from cells.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = True
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bam file', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output count table file", required=NOT_DEBUG)
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    args.input="/scratch/jbrown_lab/shengq2/projects/20190930_scRNA_3804_mouse_star_fix/star/result/LD2_Wild_1824.bam"
    args.output="/scratch/jbrown_lab/shengq2/projects/20190930_scRNA_3804_mouse_star_fix/star/result/LD2_Wild_1824.count"
  
  logger = logging.getLogger('10xCount')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  count(args.output, args.input, logger)
  
if __name__ == "__main__":
    main()
