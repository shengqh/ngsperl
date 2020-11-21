import argparse
import sys
import logging
import os
import pybedtools 
import requests
import shutil

def runCmd(cmd, logger):
  logger.info(cmd)
  os.system(cmd)

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Get active gene",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input peak bed file', required=NotDEBUG)
parser.add_argument('-g', '--genome', action='store', nargs='?', help='Input genome (mm10/hg19/hg38)', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output prefix", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/jbrown_lab/shengq2/projects/20200720_chipseq_GSE140641_mouse/croo/result/H3K27ac_activated/peak/overlap_reproducibility/overlap.optimal_peak.narrowPeak.gz.bed"
  args.genome = "mm10"
  args.output = "H3K27ac_activated"

logger = logging.getLogger('activeGenes')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

if not os.path.isfile(args.genome):
  annotation_file="%s_refseq.ucsc" % args.genome
  if not os.path.isfile(annotation_file):
    logger.info("download %s ..." % url)
    url="https://raw.githubusercontent.com/linlabcode/pipeline-tools/master/pipeline_tools/annotation/%s_refseq.ucsc" % args.genome
    download_command = "wget " + url
    runCmd(download_command, logger)
else:
  annotation_file=args.genome

region_file = '%s_promoter_regions.bed' % args.output
logger.info("convert to %s ..." % region_file)
with open(annotation_file, "rt") as infile, open(region_file, 'wt') as outfile:
  infile.readline()
  for line in infile:
    values = line.split('\t')
    chrom, strand, start, stop, name = values[2:6] + values[12:13]
    choord = int(start) if strand == '+' else int(stop)
    outfile.write('\t'.join([chrom, str(max(0, choord - 1000)), str(choord + 1000), '.', '', strand, name]) + '\n')

logger.info("generate active genes ...")

tmp_bed = args.output + ".tmp.bed"
method_command = "bedtools intersect -a \"%s\" -b \"%s\" > %s" % (region_file, args.input, tmp_bed)
runCmd(method_command, logger)

intersect = pybedtools.BedTool(tmp_bed)

# promoters = pybedtools.BedTool(tmp_bed)
# bed = pybedtools.BedTool(args.input)
# intersect = promoters.intersect(bed)

gene_list = [line[6] for line in intersect]
gene_list = sorted(list(set(gene_list)))

final_file = args.output + ".TSS_ACTIVE_-1000_1000.txt"
logger.info(f"writing result to {final_file} ...")
with open(final_file, 'wt') as outfile:
  outfile.write('\n'.join(gene_list))

#os.remove(annotation_file)
os.remove(region_file)
os.remove(tmp_bed)

logger.info("done.")