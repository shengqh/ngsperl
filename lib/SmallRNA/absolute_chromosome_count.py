import logging
import argparse
import pysam
from DupCountUtils import readDupCountMap

DEBUG = 0

if DEBUG:
  inputFile="/scratch/vickers_lab/projects/20210504_6282_RA_smRNAseq_mouse_byTiger/covid19/bowtie1_covid19/result/WTwdM19.bam"
  countFile="/scratch/vickers_lab/projects/20210504_6282_RA_smRNAseq_mouse_byTiger/preprocessing/identical/result/WTwdM19_clipped_identical.fastq.dupcount"
  outputFile="/scratch/vickers_lab/projects/20210504_6282_RA_smRNAseq_mouse_byTiger/covid19/bowtie1_covid19_count/result/WTwdM19.count"
  name_map_file=None
else:
  parser = argparse.ArgumentParser(description="Extract mapped reads to Fastq.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bam file')
  parser.add_argument('-n', '--name_map', action='store', nargs='?', help='Input chromosome to name map file')
  parser.add_argument('-c', '--count', action='store', nargs='?', help="Original dupcount file")
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output result file")

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  countFile = args.count
  outputFile = args.output
  name_map_file = args.name_map

logger = logging.getLogger('count')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("Reading counts from " + countFile)
counts = readDupCountMap(countFile)

if name_map_file == None:
  name_map = None
else:
  name_map = {}
  with open(name_map_file, "r") as sr:
    sr.readline()
    for line in sr:
      parts = line.rstrip().split('\t')
      chromosome = parts[0]
      name = int(parts[1])
      name_map[chromosome] = name


count_map = {}
logger.info("Reading BAM file " + inputFile)
with pysam.Samfile(inputFile) as sam:
  for read in sam.fetch():
    refname = read.reference_name
    count_map.setdefault(refname, {})[read.query_name] = counts[read.query_name]

with open(outputFile, "w") as sw:
  sw.write("Chrom\tCount\n")
  for chrom in sorted(count_map.keys()):
    chrom_map = count_map[chrom]
    chrom_count = sum(chrom_map.values())
    sw.write(f"{chrom}\t{chrom_count}\n")

logger.info("Done")
