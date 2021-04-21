import argparse
import sys
import logging
import os
import subprocess

def bamsnap(logger, bed_file, bam_list_file, output_file, bamsnap_option_file, no_gene_track, extend_bases):
  bamsnap_option = ""
  with open(bamsnap_option_file, "r") as fin:
    for line in fin:
      parts = line.rstrip().split('\t')
      if parts[0] != '':
        bamsnap_option = bamsnap_option + " " + parts[1] + " " + parts[0]
  logger.info("bamsnap_option=" + bamsnap_option)

  gene_locus_map = {}
  with open(bed_file, "r") as ins:
    for line in ins:
      parts = line.strip().split('\t')
      if extend_bases > 0:
        locus = "%s:%d-%d" % (parts[0], int(parts[1]) - extend_bases, int(parts[2]) + extend_bases)
      else:
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

  gene_option = "" if no_gene_track else "gene "

  gene_file_map={}
  for gene in sorted(gene_locus_map.keys()):
    locus = gene_locus_map[gene]
    logger.info("Drawing " + gene + " " + locus + " ...")
    pngfile = gene + ".png"
    snapCommand = "bamsnap %s -draw coordinates bamplot %s -bamplot coverage -width 1000 -height %d -pos %s -out %s %s %s" % (bamsnap_option, gene_option, height, locus, pngfile, bam_name_str, bam_file_str)
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
  parser.add_argument('-b', '--bam', action='store', nargs='?', help="Bam list file", required=NOT_DEBUG)
  parser.add_argument('-c', '--bamsnap_option_file', action='store', nargs='?', help="bamsnap option file")
  parser.add_argument('-e', '--extend_bases', action='store', type=int, default=1000, nargs='?', help="Extending X bases before and after coordinates")
  parser.add_argument('--no_gene_track', action='store_true', help="No gene track")
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output list file", required=NOT_DEBUG)
  
  args = parser.parse_args()
  
  if DEBUG:
    args.input="/scratch/jbrown_lab/shengq2/projects/20210413_cutrun_6133_human/annotation_genes_locus/result/cutrun_6133.shift10000.bed"
    args.bam="/scratch/jbrown_lab/shengq2/projects/20210413_cutrun_6133_human/annotation_genes_bamsnap/result/cutrun_6133__fileList1.list"
    args.bamsnap_option_file="/scratch/jbrown_lab/shengq2/projects/20210413_cutrun_6133_human/annotation_genes_bamsnap/result/cutrun_6133__fileList2.list"
    args.no_gene_track=False
    args.output="/scratch/jbrown_lab/shengq2/projects/20210413_cutrun_6133_human/annotation_genes_bamsnap/result/cutrun_6133.txt"
  
  logger = logging.getLogger('bamsnap')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  bamsnap(logger, args.input, args.bam, args.output, args.bamsnap_option_file, args.no_gene_track, args.extend_bases)
  
if __name__ == "__main__":
    main()
