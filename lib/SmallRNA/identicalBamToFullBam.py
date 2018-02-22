import pysam
import sys
import gzip
import os
import logging
import argparse
from DupCountUtils import readDupCountMap 

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  if DEBUG:
    inputFile="/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_genome/bowtie1_bacteria_group2_pm/result/Urine_WT_14/Urine_WT_14.bam"
    countFile = "/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_unmapped_reads/result/Urine_WT_14_clipped_identical.unmapped.fastq.dupcount"
    chromosome = "NC_012660.1"
    outputFile="/scratch/cqs/shengq2/temp/Urine_WT_14.bam"
  else:
    parser = argparse.ArgumentParser(description="Generate full BAM from identical BAM.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
    parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bam file')
    parser.add_argument('-c', '--count', action='store', nargs='?', help="Input dupcount file")
    parser.add_argument('-r', '--chromosome', action='store', nargs='?', help="Filter chromosome")
    parser.add_argument('-o', '--output', action='store', nargs='?', help="Output bam file")
  
    args = parser.parse_args()
    
    print(args)
    
    inputFile = args.input
    countFile = args.count
    chromosome = args.chromosome
    outputFile = args.output
  
  logger = logging.getLogger('bamToBam')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  logger.info("Reading count file " + countFile + " ...")  
  readMap = readDupCountMap(countFile)
  
  unsorted = outputFile + ".unsorted.bam"
  with pysam.Samfile(inputFile) as sam:
    with pysam.AlignmentFile(unsorted, "wb", template=sam) as outf:
      for read in sam.fetch(chromosome):
        query_name = read.query_name
        query_count = readMap[query_name]
        for idx in range(1, query_count):
          read.query_name = query_name + ":" + str(idx)
          outf.write(read)
     
  
  filename, file_extension = os.path.splitext(outputFile)
  pysam.sort("-o", outputFile, unsorted)
  pysam.index(outputFile)
  os.remove(unsorted)
  
if __name__ == "__main__":
    main()
