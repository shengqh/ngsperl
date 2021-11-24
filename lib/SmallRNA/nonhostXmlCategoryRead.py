import sys
import os
import logging
import argparse

DEBUG = False

if DEBUG:
  inputFile="/scratch/vickers_lab/projects/20210901_smallRNA_CACvsControl_DC_2845_2792_humanv5_byMars/nonhost_library/bowtie1_tRNA_pm_table/result/tRNA_pm_DC_2845_2792_human.category.count.xml"
  seqName="Bacteria"
  outputFile="/scratch/cqs/shengq2/temp/tRNA_pm_DC_2845_2792_human.Bacteria.count"
else:
  parser = argparse.ArgumentParser(description="Get read count from non host database",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count or xml file')
  parser.add_argument('-c', '--category', action='store', nargs='?', help='Input category name')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output count file")

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  seqName = args.category
  outputFile = args.output

logger = logging.getLogger('xmlCount')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

if inputFile.endswith(".count"):
  inputFile = inputFile + ".xml"

seq_count = {}
samples = set()
logger.info("Processing " + inputFile + " ...")
with open(inputFile, 'rt') as xml_file:
  for line in xml_file:
    if '</queries>' in line:
      break

    if '<query ' in line:
      for nextline in xml_file:
        if seqName in nextline:
          parts = line.split("\"")
          query_sequence = parts[3]
          query_sample = parts[5]
          query_count = int(parts[7])
          seq_count.setdefault(query_sequence, {})[query_sample] = query_count
          samples.add(query_sample)
          break
        if '</query>' in nextline:
          break

seqc=[[k, sum(seq_count[k].values())] for k in seq_count.keys()]
seqc.sort(key=lambda x:x[1], reverse=True)
sequences=[x[0] for x in seqc]

logger.info("Writing result to " + outputFile + " ...")
samples = sorted([s for s in samples])
with open(outputFile, "wt") as fout:
  fout.write("Sequence\t%s\n" % "\t".join(samples))
  for seq in sequences:
    fout.write(seq)
    count_map = seq_count[seq]
    for sample in samples:
      if sample in count_map:
        fout.write(f"\t{count_map[sample]}")
      else:
        fout.write(f"\t0")
    fout.write('\n')


logger.info("done")