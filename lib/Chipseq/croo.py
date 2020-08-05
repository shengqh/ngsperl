import argparse
import sys
import logging
import os
import csv
import gzip
  
DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Perform croo to retrive wdl result.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input wdl result folder', required=NotDEBUG)
parser.add_argument('-n', '--name', action='store', nargs='?', help='Input sample name', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output folder", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/jbrown_lab/shengq2/projects/20200720_chipseq_GSE140641_mouse/Encode/result/H3K27ac_activated/chip"
  args.name = "H3K27ac_activated"
  args.output = "/scratch/jbrown_lab/shengq2/projects/20200720_chipseq_GSE140641_mouse/croo/result/H3K27ac_activated"

logger = logging.getLogger('croo')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

subfolders = [ f.path for f in os.scandir(args.input) if f.is_dir() ]
metafiles = [os.path.join(sf, "metadata.json") for sf in subfolders if os.path.exists(os.path.join(sf, "metadata.json"))]

if len(metafiles) > 1:
  raise Exception("Multiple metadata.json found: %s" % ",".join(metafiles))
elif len(metafiles) == 0:
  raise Exception("No metadata.json found: %s" % args.input)

cmd = "croo --out-dir %s %s" % (args.output, metafiles[0])
logger.info(cmd)
os.system(cmd)

rep1bam = "%s/align/rep1/%s.nodup.bam" % (args.output, args.name)
ctl1bam = "%s/align/ctl1/input.nodup.bam" % (args.output)
peak = "%s/peak/overlap_reproducibility/overlap.optimal_peak.narrowPeak.gz" % args.output

cmd = "samtools index %s" % rep1bam
logger.info(cmd)
os.system(cmd)

cmd = "samtools index %s" % ctl1bam
logger.info(cmd)
os.system(cmd)

cmd = "gunzip -c %s > %s.bed" % (peak, peak)
logger.info(cmd)
os.system(cmd)

logger.info("done.")

