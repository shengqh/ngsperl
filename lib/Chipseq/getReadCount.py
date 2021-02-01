import argparse
import sys
import logging
import os
import shutil
import math
from collections import OrderedDict
import pysam
import csv

class LocusItem(object):
  def __init__(self, chromosome, TSS_start, TSS_end, name, score, strand, gene_body_start, gene_body_end):
    self.Chromosome = chromosome
    self.TSS_start = TSS_start
    self.TSS_end = TSS_end
    self.Name = name
    self.Score = score
    self.Strand = strand
    self.Gene_body_start = gene_body_start
    self.Gene_body_end = gene_body_end
    self.Locus = "%s:%d-%d" % (chromosome, TSS_start, TSS_end)
    self.Count = {}

def readGenes(deseq2sigFile):
  result = {}
  with open(deseq2sigFile, newline='') as csvfile:
    spamreader = csv.DictReader(csvfile)
    for row in spamreader:
      result[row["Feature_gene_name"]] = float(row["log2FoldChange"]) > 0
  return(result)

def readBedFile(fileName):
  result = []
  with open(fileName, "r") as fin:
    for line in fin:
      parts = line.rstrip().split('\t')
      item = LocusItem(parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4], parts[5], int(parts[6]), int(parts[7]))
      result.append(item)
  return(result)

def getCount(fsam, chrom, start, end):
  query_names = set()
  for read in fsam.fetch(chrom, start, end):
    query_names.add(read.query_name)
  return(len(query_names))

def runCmd(cmd, logger):
  logger.info(cmd)
  os.system(cmd)

def getRPM(count, totalCount):
  return(count * 1000000.0 / totalCount)

def getRatio(sample_rpm, control_rpm):
  if control_rpm == 0:
    return(50)
  elif (sample_rpm == 0):
    return(-50)
  else:
    return(math.log2(sample_rpm / control_rpm))

def getReadCount(logger, output_file, gene_file, bed_file, sample_bam, control_bam):
  if os.path.isfile(output_file):
    return

  genes = readGenes(gene_file)
  logger.info("Total %d genes from gene file" % len(genes))

  all_items = readBedFile(bed_file)
  items = [item for item in all_items if item.Score != 0 and item.Name in genes]
  logger.info("Total %d valid genes from bed file" % len(items))

  logger.info("reading " + sample_bam + " ...")
  with open(output_file, "wt") as fout:
    with pysam.Samfile(sample_bam) as fin, pysam.Samfile(control_bam) as fcon:
      scount = fin.count()
      ccount = fcon.count()
      with open(output_file + ".totalcount", "wt") as ftout:
        ftout.write("Sample\t%d\nControl\t%d\n" % (scount, ccount) )

      logger.info("Sample read count=%d, control read count=%d" % (scount, ccount))
      logger.info("start counting %s ..." % bed_file)
      fout.write("Gene\tUpregulated\tSampleTSSCount\tSampleTSSRPM\tSampleGenebodyCount\tSampleGenebodyRPM\tControlTSSCount\tControlTSSRPM\tControlGenebodyCount\tControlGenebodyRPM\tlog2TSSRatio\tlog2GenebodyRatio\n")
      for item in items:
        s_tss = getCount(fin, item.Chromosome, item.TSS_start, item.TSS_end)
        s_gb = getCount(fin, item.Chromosome, item.Gene_body_start, item.Gene_body_end)
        c_tss = getCount(fcon, item.Chromosome, item.TSS_start, item.TSS_end)
        c_gb = getCount(fcon, item.Chromosome, item.Gene_body_start, item.Gene_body_end)
        s_tss_rpm = getRPM(s_tss, scount)
        s_gb_rpm = getRPM(s_gb, scount)
        c_tss_rpm = getRPM(c_tss, ccount)
        c_gb_rpm = getRPM(c_gb, ccount)
        tss_ratio = getRatio(s_tss_rpm, c_tss_rpm)
        gb_ratio = getRatio(s_gb_rpm, c_gb_rpm)

        fout.write("%s\t%d\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%f\t%f\n" %
          (item.Name,
          genes[item.Name],
          s_tss,
          s_tss_rpm,
          s_gb,
          s_gb_rpm,
          c_tss,
          c_tss_rpm,
          c_gb,
          c_gb_rpm,
          tss_ratio,
          gb_ratio))
        logger.info("%s: tss_ratio=%f, gb_ratio=%f" % (item.Name, tss_ratio, gb_ratio))

    logger.info("done.")

def getLogger():
  logger = logging.getLogger('getReadCount')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  return(logger)

def main():
  DEBUG=True
  NotDEBUG=not DEBUG

  parser = argparse.ArgumentParser(description="Get read count",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bam file', required=NotDEBUG)
  parser.add_argument('-g', '--gene_file', action='store', nargs='?', help='Input gene file', required=NotDEBUG)
  parser.add_argument('-b', '--bed_file', action='store', nargs='?', help='Input coordinate bed file', required=NotDEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output prefix", required=NotDEBUG)

  args = parser.parse_args()
  if DEBUG:
    #args.input = "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/TNF_Veh_886.noChrM.bam"
    #args.output = "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/TNF_Veh_886.tss_genebody.tsv"
    #args.control = "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/No_Treatment_886.noChrM.bam"
    args.input = "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/THZ_TNF_886.noChrM.bam"
    args.output = "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/THZ_TNF_886_vs_TNF_Veh.tss_genebody.tsv"
    args.control = "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/TNF_Veh_886.noChrM.bam"
    args.gene_file = "/scratch/jbrown_lab/shengq2/projects/20201208_rnaseq_464_hg38_endothelial_cell/deseq2_proteincoding_genetable/result/TNF_THZ25_vs_TNF_Veh_min5_fdr0.05_DESeq2_sig_genename.txt"
    args.bed_file = "/scratch/cqs_share/references/gencode/GRCh38.p13/gencode.v36.annotation.tss_genebody.bed"

  logger = getLogger()
  getReadCount(logger, args.output, args.gene_file, args.bed_file, args.input, args.control)

if __name__ == "__main__":
  #main()
  logger = getLogger()

  deseq2file="/scratch/jbrown_lab/shengq2/projects/20201208_rnaseq_464_hg38_endothelial_cell/deseq2_proteincoding_genetable/result/TNF_Veh_vs_No_Treatment_min5_fdr0.05_DESeq2_sig.csv"

  # getReadCount(logger, 
  #         "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/THZ_TNF_886_vs_TNF_Veh.tss_genebody.tsv", 
  #         deseq2file,
  #         "/scratch/cqs_share/references/gencode/GRCh38.p13/gencode.v36.annotation.tss_genebody.bed", 
  #         "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/THZ_TNF_886.noChrM.bam", 
  #         "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/TNF_Veh_886.noChrM.bam")
  
  getReadCount(logger, 
          "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/TNF_Veh_886_vs_No_Treatment.tss_genebody.tsv", 
          deseq2file,
          "/scratch/cqs_share/references/gencode/GRCh38.p13/gencode.v36.annotation.tss_genebody.bed", 
          "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/TNF_Veh_886.noChrM.bam", 
          "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/No_Treatment_886.noChrM.bam")
  
  getReadCount(logger, 
          "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/THZ_TNF_886_vs_No_Treatment.tss_genebody.tsv", 
          deseq2file,
          "/scratch/cqs_share/references/gencode/GRCh38.p13/gencode.v36.annotation.tss_genebody.bed", 
          "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/THZ_TNF_886.noChrM.bam", 
          "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/No_Treatment_886.noChrM.bam")
  
  
  getReadCount(logger, 
          "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/TNF_Veh_485_vs_No_Treatment.tss_genebody.tsv", 
          deseq2file,
          "/scratch/cqs_share/references/gencode/GRCh38.p13/gencode.v36.annotation.tss_genebody.bed", 
          "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/TNF_Veh_485.noChrM.bam", 
          "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/No_Treatment_485.noChrM.bam")
  
  getReadCount(logger, 
          "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/THZ_TNF_485_vs_No_Treatment.tss_genebody.tsv", 
          deseq2file,
          "/scratch/cqs_share/references/gencode/GRCh38.p13/gencode.v36.annotation.tss_genebody.bed", 
          "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/THZ25_TNF_485.noChrM.bam", 
          "/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/No_Treatment_485.noChrM.bam")
