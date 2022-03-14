import logging
import re
import argparse
import gzip

def initialize_logger():
  logger = logging.getLogger('count')
  loglevel = logging.INFO
  logger.setLevel(loglevel)

  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')    
 
  # create console handler and set level to info
  handler = logging.StreamHandler()
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  return(logger)

logger = initialize_logger()

DEBUG = False
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="Get read count in chromosomes",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file', required=NOT_DEBUG)
parser.add_argument('-s', '--species', action='store', nargs='?', help='Input species file', required=NOT_DEBUG)
parser.add_argument('-c', '--count', action='store', nargs='?', help='Input count file', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output summary file", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="/scratch/vickers_lab/projects/20220301_4893_2_RA_smRNA_mouse_v5_bacteria/bowtie_count/result/bacteria__fileList1.list"
  args.species="/data/cqs/references/bacteria/20220228_bacteria.species.txt"
  args.count='/scratch/vickers_lab/projects/20200625_4893_2_RA_smRNA_mouse_v5_byTiger/intermediate_data/bowtie1_genome_unmapped_reads/result/WT_SH_5_clipped_identical.unmapped.fastq.dupcount'
  args.output="/scratch/vickers_lab/projects/20220301_4893_2_RA_smRNA_mouse_v5_bacteria/bowtie_count/result/bacteria.count"

count_map={}
logger.info(f"reading count file {args.count}")
with open(args.count, "rt") as fin:
  fin.readline()
  for line in fin:
    parts=line.split('\t')
    count_map[parts[0]] = parts[1]

logger.info(f"reading species file {args.species}")
species_map = {}
with open(args.species, "rt") as fin:
  fin.readline()
  for line in fin:
    parts = line.rstrip().split('\t')
    species_map[parts[0]] = parts[1]

read_map = {}
with open(args.input, "rt") as fl:
  for line in fl:
    parts = re.split('\s+', line.rstrip())
    bowtie_file = parts[0]
    logger.info(f"parsing {bowtie_file}")
    with gzip.open(bowtie_file, "rt") as fin:
      for bl in fin:
        bparts = bl.split('\t')
        query = bparts[0].split(' ')[0]
        species = species_map[bparts[2]]
        read_map.setdefault(query,set()).add(species)

logger.info(f"merge all bowtie result ...")
all_queries = [[query, int(count_map[query]), ",".join(read_map[query])] for query in read_map.keys()]   
all_queries.sort(key=lambda x:x[1], reverse=True)

logger.info(f"output to {args.output} ...")
with gzip.open(args.output, "wt") as fout:
  fout.write("read\tcount\tspecies\n")
  for query in all_queries:
    fout.write(f"{query[0]}\t{query[1]}\t{query[2]}\n")

logger.info("done")
