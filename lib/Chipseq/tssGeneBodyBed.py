import argparse
import sys
import logging
import os
import shutil
from collections import OrderedDict
from gtfparse import read_gtf

class Gene:
  def __init__(self, chrom, start, end, name, strand):
    self.chrom = chrom
    self.start = start
    self.end = end
    self.name = name
    self.strand = strand
    self.overlap = False
    self.calc_tss_genebody()

  def calc_tss_genebody(self):
    if self.strand == '+':
      self.tss_start = self.start - 300
      self.tss_end = self.start + 300
      self.gb_start = self.start + 301
      self.gb_end = self.end + 3500
      self.all_start = self.tss_start
      self.all_end = self.gb_end
    else:
      self.tss_start = self.end - 300
      self.tss_end = self.end + 300
      self.gb_start = max(1, self.start - 3500)
      self.gb_end = self.end - 301
      self.all_start = self.gb_start
      self.all_end = self.tss_end

  def output_tss(self, fout):
    fout.write("%s\t%d\t%d\t%s\t%d\t%s\n" % (self.chrom, self.tss_start, self.tss_end, self.name, 0 if self.overlap else 1000, self.strand))

  def output_gene_body(self, fout):
    fout.write("%s\t%d\t%d\t%s\t%d\t%s\n" %(self.chrom, self.gb_start, self.gb_end, self.name, 0 if self.overlap else 1000, self.strand))

  def output_tss_gene_body(self, fout):
    fout.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\n" %(self.chrom, self.tss_start, self.tss_end, self.name, 0 if self.overlap else 1000, self.strand, self.gb_start, self.gb_end))

  def mark_overlap(self, another):
    if self.all_end < another.all_start:
      return

    if another.all_end < self.all_start:
      return

    self.overlap = True
    another.overlap = True

def runCmd(cmd, logger):
  logger.info(cmd)
  os.system(cmd)

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Get bed file for TSS and gene body",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input gtf file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output tss gene bed file", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/cqs_share/references/gencode/GRCh38.p13/gencode.v36.annotation.gtf"
  args.output = "/scratch/cqs_share/references/gencode/GRCh38.p13/gencode.v36.annotation"

logger = logging.getLogger('tssGeneBody')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("reading " + args.input + " ...")
chrom_genes = OrderedDict()
with open(args.input, "rt") as fin:
  for line in fin:
    if line.startswith("#"):
      continue

    parts = line.split('\t')
    if parts[2] != 'gene':
      continue
    
    if not 'gene_type "protein_coding"' in parts[8]:
      continue

    infos = parts[8].split('; ')
    for info in infos:
      if info.startswith("gene_name"):
        gene_name = info[11:][:-1]
        start = int(parts[3])
        end = int(parts[4])
        gene = Gene(parts[0], start, end, gene_name, parts[6])
        chrom_genes.setdefault(gene.chrom, []).append(gene)

with open(args.output + ".tss.bed", "wt") as ftss:
  with open(args.output + ".genebody.bed", "wt") as fbody:
    with open(args.output + ".tss_genebody.bed", "wt") as fall:
      for chrom in chrom_genes.keys():
        cur_genes = chrom_genes[chrom]
        logger.info("checking overlap in chrom %s ..." % chrom)
        for i in range(0, len(cur_genes)):
          for j in range(i+1, len(cur_genes)):
            cur_genes[i].mark_overlap(cur_genes[j])

        for cur_gene in cur_genes:
          cur_gene.output_tss(ftss)
          cur_gene.output_gene_body(fbody)
          cur_gene.output_tss_gene_body(fall)

logger.info("done.")