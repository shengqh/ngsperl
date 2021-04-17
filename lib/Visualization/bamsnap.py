import argparse
import sys
import logging
import os
import subprocess

def bamsnap(logger, bed_file, bam_list_file, output_file, discard_gene=True, refversion="hg38"):
  gene_locus_map = {}
  with open(bed_file, "r") as ins:
    for line in ins:
      parts = line.strip().split('\t')
      locus = "%s:%s-%s" % (parts[0], parts[1], parts[2])
      gene = parts[4]
      gene_locus_map[gene] = locus

  bam_file_map = {}
  with open(bam_list_file, "r") as ins:
    for line in ins:
      parts = line.strip().split('\t')
      bam_file_map[parts[1]] = parts[0]

  bam_names = sorted(bam_file_map.keys())
  bam_name_str = '-title ' + '"' + '" "'.join(bam_names) + '"'
  bam_file_str = '-bam ' + ' '.join([bam_file_map[n] for n in bam_names])

  height = max(len(bam_names) * 300 + 1100, 2000)

  draw_gene_option = "" if discard_gene else "gene "
  gene_file_map={}
  for gene in sorted(gene_locus_map.keys()):
    locus = gene_locus_map[gene]
    logger.info("Drawing " + gene + " " + locus + " ...")
    pngfile = gene + ".png"
    snapCommand = "bamsnap -draw coordinates bamplot " + draw_gene_option + " -bamplot coverage -refversion %s -width 1000 -height %d -pos %s -out %s %s %s" % (refversion, height, locus, pngfile, bam_name_str, bam_file_str)
    print(snapCommand)
    subprocess.run(snapCommand, check=True, shell=True)
    gene_file_map[gene] = pngfile

  with open(output_file, "wt") as fout:
    for gene in sorted(gene_file_map.keys()):
      fout.write(gene_file_map[gene] + "\t" + gene + "\n")

  logger.info("done.")
        
def main():
  parser = argparse.ArgumentParser(description="Call bamsnap for gene visualization",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bed file', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output list file", required=NOT_DEBUG)
  parser.add_argument('-b', '--bam', action='store', nargs='?', help="Bam list file", required=NOT_DEBUG)
  parser.add_argument('-d', '--discard_gene', action='store_true', help="Discard gene for genome not hg38/hg19/mm10")
  parser.add_argument('--refversion', action='store', nargs='?', help="Genome version (hg38/hg19/mm10)", required=NOT_DEBUG)
  parser.add_argument('--coverage_color', action='store', nargs='?', help="Hex color for coverage (000000 for black)")
  
  args = parser.parse_args()
  
  if DEBUG:
    args.input="/scratch/jbrown_lab/shengq2/projects/20210321_cutrun_6048_human/annotation_genes_locus/result/cutrun_6048.bed"
    args.bam="/scratch/jbrown_lab/shengq2/projects/20210321_cutrun_6048_human/annotation_genes_bamsnap/result/cutrun_6048__fileList1.list"
    args.output="/scratch/jbrown_lab/shengq2/projects/20210321_cutrun_6048_human/annotation_genes_bamsnap/result/cutrun_6048.csv"
    args.discard_gene=False
    args.refversion="hg38"
  
  logger = logging.getLogger('bamsnap')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  bamsnap(logger, args.input, args.bam, args.output, args.discard_gene, args.refversion)
  
if __name__ == "__main__":
    main()
