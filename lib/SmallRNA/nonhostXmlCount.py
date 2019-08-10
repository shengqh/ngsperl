import sys
import os
import logging
import argparse

DEBUG = False

if DEBUG:
  listFile="/data/stein_lab/mjo_sRNA_data/20170206_michelle_smallRNA_2868_human_nonhost_genome_count/result/2868__fileList1.list"
  outputFile="/data/stein_lab/mjo_sRNA_data/20170206_michelle_smallRNA_2868_human_nonhost_genome_count/result/2868.microbial.tsv"
else:
  parser = argparse.ArgumentParser(description="Get unique read count from overlapped count result.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count xml file list')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output fastq file")

  args = parser.parse_args()
  
  print(args)
  
  listFile = args.input
  outputFile = args.output

logger = logging.getLogger('xmlCount')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

inputFiles = []
with open(listFile, 'r') as sr:
  for line in sr:
    parts = line.rstrip().split('\t')
    inputFiles.append(parts[0])

sampleCount = {}
mapped = set()
for inputFile in inputFiles:
  logger.info("Processing " + inputFile + " ...")
  with open(inputFile, 'rt') as xml_file:
    for line in xml_file:
      if '</queries>' in line:
        break

      if '<query ' in line:
        parts = line.split("\"")
        query_name = parts[1]
        if query_name not in mapped:
          query_count = int(parts[7])
          query_sample = parts[5]
          mapped.add(query_name)

          if query_sample not in sampleCount:
            sampleCount[query_sample] = query_count
          else:
            sampleCount[query_sample] = sampleCount[query_sample]  + query_count

mapped.clear()

logger.info("Writing result to " + outputFile + " ...")
samples = sorted(sampleCount.keys())
with open(outputFile, "wt") as fout:
  fout.write("Sample\tCount\n")
  for sample in samples:
    fout.write("%s\t%s\n" % (sample, sampleCount[sample]))
  
