import logging
import argparse
import pysam
import pickle
from DupCountUtils import readDupCountQueries

DEBUG=False

if DEBUG:
  inputFile="/scratch/vickers_lab/projects/20210817_smallRNA_contamination/3626_RA/bowtie1_contamination_all_01_align/result/HDL_36_JPN.bam"
  countFile="/scratch/vickers_lab/projects/20210818_3626_RA_smRNA_human_contamination/preprocessing/identical/result/HDL_36_JPN_clipped_identical.fastq.dupcount"
  outputFile="/scratch/vickers_lab/projects/20210817_smallRNA_contamination/3626_RA/bowtie1_contamination_all_02_count/result/HDL_36_JPN.count"
  name_map_file=None
  extend_bam=False
else:
  parser = argparse.ArgumentParser(description="Extract mapped reads to Fastq.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bam file')
  parser.add_argument('-n', '--name_map', action='store', nargs='?', help='Input chromosome to name map file')
  parser.add_argument('-c', '--count', action='store', nargs='?', help="Original dupcount file")
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output result file")
  parser.add_argument('-e', '--extend_bam', action='store_true', help="Output extend bam file")

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  countFile = args.count
  outputFile = args.output
  name_map_file = args.name_map
  extend_bam = args.extend_bam

logger = logging.getLogger('count')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("Reading counts from " + countFile)
queries = readDupCountQueries(countFile, 0)
query_map = {q.Name:q for q in queries}

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
    query = query_map[read.query_name]
    count_map.setdefault(refname, {})[query.Sequence] = [query.Count, read.reference_start, read.reference_end]

with open(outputFile, "w") as sw:
  sw.write("Chrom\tUniqueCount\tReadCount\n")
  for chrom in sorted(count_map.keys()):
    chrom_map = count_map[chrom]
    unique_count = len(chrom_map)
    chrom_count = sum([v[0] for v in chrom_map.values()])
    sw.write(f"{chrom}\t{unique_count}\t{chrom_count}\n")
  sw.write(f"total\t{len(query_map)}\t{sum([q.Count for q in query_map.values()])}\n")

with open(outputFile + ".obj", 'wb') as fout:
  pickle.dump(count_map, fout)  

chroms = sorted(count_map.keys())
seqs={}
for chrom in chroms:
  chrom_map = count_map[chrom]
  items = sorted(chrom_map.items(), key=lambda item: item[1][0],  reverse=True)
  for idx in range(0, min(20, len(items))):
    item=items[idx]
    seqs[item[0]] = item[1][0]

seq_items = sorted(seqs.items(), key=lambda item: item[1],  reverse=True)
#print(seq_items)

with open(outputFile + ".seq", "w") as sw:
  sw.write("Index\tSequence\tCount\t%s\n" % "\t".join(chroms))
  for idx in range(0, len(seq_items)):
    item = seq_items[idx]
    seq = item[0]
    count = item[1]
    sw.write(f"{idx+1}\t{seq}\t{count}")
    for chrom in chroms:
      chrom_map = count_map[chrom]
      if seq in chrom_map:
        citem = chrom_map[seq]
        sw.write(f"\t{citem[1]}-{citem[2]}")
      else:
        sw.write("\t")
    sw.write("\n")

if extend_bam:
  extend_bam_file = outputFile + ".bam"
  with pysam.Samfile(inputFile) as sam:
    header = sam.header
    with pysam.AlignmentFile(extend_bam_file, "wb", header=header) as outf:
      for read in sam.fetch():
        query_name = read.query_name
        query_count = query_map[query_name].Count
        for idx in range(0, query_count):
          read.query_name = query_name + ":" + str(idx)
          outf.write(read)
    
  pysam.index(extend_bam_file)

logger.info("Done")
