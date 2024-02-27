import sys
import logging
import os
import argparse
import pysam
import numpy as np
import pickle 
from collections import OrderedDict
from asyncore import read

class Item(object):
  def __init__(self, chrom, start, end, category, name, strand):
    self.chrom = chrom
    self.start = start
    self.end = end
    self.category = category
    self.name = name
    self.strand = strand
  
  def __str__(self):
    return f"{self.chrom}\t{self.start}\t{self.end}\t{self.name}\t{self.category}\t{self.strand}"
  
def gtf_to_items(gtf, logger):
  logger.info("Processing " + gtf + " ...")

  items = []
  with open(gtf, "rt") as fin:
    count = 0
    lastitem = None
    for line in fin:
      if (line.startswith('#')):
        continue
      parts = line.rstrip().split('\t')
      #print(parts)
      if parts[2] == 'transcript':
        gene_chrom = parts[0]
        gene_start = int(parts[3]) - 1
        gene_end = int(parts[4])
        strand = parts[6]
        transcript_id = parts[8].split('transcript_id "', 1)[1].split('"', 1)[0]
        curitem = Item(gene_chrom, gene_start, gene_end, "transcript", transcript_id, strand)
        #print(f"transcript={curitem}")
        count += 1
        lastitem = curitem
      elif parts[2] == 'exon':
        gene_chrom = parts[0]
        gene_start = int(parts[3]) - 1
        gene_end = int(parts[4])
        transcript_id = parts[8].split('transcript_id "', 1)[1].split('"', 1)[0]
        strand = parts[6]
        curitem = Item(gene_chrom, gene_start, gene_end, "exon", transcript_id, strand)
        if not lastitem is None:
          if lastitem.category == 'exon' and lastitem.name == curitem.name:
            if curitem.strand == '+':
              iitem = Item(curitem.chrom, lastitem.end, curitem.start, "intron", curitem.name, curitem.strand)
            else:
              iitem = Item(curitem.chrom, curitem.end, lastitem.start, "intron", curitem.name, curitem.strand)
            items.append(iitem)
          #print(f"{iitem}")
        items.append(curitem)
        lastitem = curitem
        #print(f"{curitem}")

  chrom_map = {}
  index = 0
  for item in items:
    if not item.chrom in chrom_map:
      index += 1
      chrom_map[item.chrom] = index

  items = sorted(items, key = lambda x: (chrom_map[x.chrom], x.start))

  return(items)

def merge_bed(items, logger):
  i = 0
  while i < len(items):
    ii = items[i]
    j = i + 1
    while j < len(items):
      ij = items[j]
      if ii.chrom != ij.chrom or ii.end < ij.start:
        break
      ii.end = max(ii.end, ij.end)
      #print(f"{ii} + {ij}")
      del items[j]
      #print(f"  = {ii}")
    i += 1
    if i % 10000 == 0:
      logger.info(i)

def get_chrom_len(bamfile):
  result = {}
  with pysam.AlignmentFile(bamfile, "rb") as sf:
    references = sf.references
    lengths = [ sf.get_reference_length(r) for r in sf.references ]
    for rl in zip(references, lengths):
      result[rl[0]] = rl[1]
  return(result)

def build_exon_intron_map(chrom_len_map, exons, introns, filename, logger):
  logger.info("build exon intron map ...")

  chrom_map = {}
  for item in exons:
    chrom_map[item.chrom] = chrom_len_map[item.chrom]
  for item in introns:
    chrom_map[item.chrom] = chrom_len_map[item.chrom]

  result = {}
  for chrom in chrom_map.keys():
    logger.info(f" init {chrom}")
    result[chrom] = np.array([0] * chrom_map[chrom])
  
  logger.info(f" process intron")
  count = 0
  for intron in introns:
    count += 1
    if count % 10000 == 0:
      logger.info(count)
    chrom_array = result[intron.chrom]
    chrom_array[intron.start:intron.end] = [1] * (intron.end - intron.start)

  logger.info(f" process exon")
  count = 0
  for exon in exons:
    count += 1
    if count % 10000 == 0:
      logger.info(count)
    chrom_array = result[exon.chrom]
    chrom_array[exon.start:exon.end] = [2] * (exon.end - exon.start)

  with open(filename, "wb") as fout:
    pickle.dump(result, fout, pickle.HIGHEST_PROTOCOL)

  return(result)

def main():
  parser = argparse.ArgumentParser(description="Quality coontrol of RNAseq BAM file.",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  DEBUG=False
  NOT_DEBUG=not DEBUG

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file (use "-" as stdin)', required=NOT_DEBUG)
  parser.add_argument('-g', '--gtf', action='store', nargs='?', help="Input gtf file", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output file name", required=NOT_DEBUG)

  args = parser.parse_args()

  if DEBUG:
    args.input = "/scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/star_featurecount/result/X6129CP21_Aligned.sortedByCoord.out.bam"
    args.gtf = "/data/cqs/references/gencode/GRCh38.p13/gencode.v38.annotation.gtf"
    args.output = "/scratch/weissvl/shengq2/20200831_human_wntpathway_rnaseq.v2/X6129CP21.chromosome.stats.txt"

  logger = logging.getLogger('rnaseqBamQC')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  chrom_len_map = get_chrom_len(args.input)
  #print(chrom_len_map)
  
  datafile = args.gtf + ".pickle"
  if os.path.exists(datafile):
    logger.info("read exon intron map ...")
    with open(datafile, "rb") as fin:
      site_map = pickle.load(fin)
  else:
    items = gtf_to_items(args.gtf, logger)

    logger.info("merge exons ...")
    exons = [item for item in items if item.category == "exon" ]
    merge_bed(exons, logger)

    logger.info("merge introns ...")
    introns = [item for item in items if item.category == "intron" ]
    merge_bed(introns, logger)

    site_map = build_exon_intron_map(chrom_len_map, exons, introns, datafile, logger)

  logger.info("process bam file ...")
  chrom_count_map = OrderedDict()
  logger.info(f"reading bam file {args.input} ..." )
  count = 0
  with pysam.AlignmentFile(args.input, "rb") as sf:
    for s in sf.fetch():
      count = count + 1
      #if count % 10000 == 0:
      #  logger.info(count)

      if count % 1000000 == 0:
        logger.info(count)
      #  break

      if s.is_unmapped:
        continue
      
      if s.reference_name not in site_map:
        break

      count_array = site_map[s.reference_name]
      v1 = count_array[s.reference_start]
      v2 = count_array[s.reference_end-1]
      v = max(v1, v2)

      varray = chrom_count_map.setdefault(s.reference_name, [0] * 3)
      varray[v] += 1

  logger.info("Writing final result to " + args.output + "...")
  with open(args.output, "wt") as fout:
    fout.write("Chromosome\tCategory\tCount\n")
    for chrom in chrom_count_map.keys():
      varray = chrom_count_map[chrom]
      fout.write(f"{chrom}\tintergenic\t{varray[0]}\n")
      fout.write(f"{chrom}\tintron\t{varray[1]}\n")
      fout.write(f"{chrom}\texon\t{varray[2]}\n")

  logger.info("done")

if __name__ == "__main__":
    main()
