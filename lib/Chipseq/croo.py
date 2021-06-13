import argparse
import glob
import logging
import os
  
DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Perform croo to retrive wdl result.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input wdl result folder', required=NotDEBUG)
parser.add_argument('-n', '--name', action='store', nargs='?', help='Input sample name', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output folder", required=NotDEBUG)
parser.add_argument('--croo', action='store', nargs='?', default="croo", help='Input croo command', required=NotDEBUG)
parser.add_argument('--out_def_json', action='store', nargs='?', default="croo", help='Input output definition JSON file for a WDL file', required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/workspace/shengq2/20210522_atacseq_6314_human_encode/encode_atacseq/result/IFNg_Rep_1/atac"
  args.name = "IFNg_Rep_1"
  args.output = "/workspace/shengq2/20210522_atacseq_6314_human_encode/encode_atacseq_croo/result/IFNg_Rep_1"
  args.croo = "singularity exec -c -B /gpfs52/data:/data,/workspace -e /data/cqs/softwares/singularity/cqs_encode.sif croo"
  args.out_def_json = "/data/cqs/softwares/encode/atac-seq-pipeline/atac.croo.v5.json"

logger = logging.getLogger('croo')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

subfolders = [ f.path for f in os.scandir(args.input) if f.is_dir() ]
metafiles = [os.path.join(sf, "metadata.json") for sf in subfolders if os.path.exists(os.path.join(sf, "metadata.json"))]

if len(metafiles) > 1:
  raise Exception("Multiple metadata.json found: %s" % ",".join(metafiles))
elif len(metafiles) == 0:
  raise Exception("No metadata.json found: %s" % args.input)

cmd = "%s --method copy --out-dir %s %s" % (args.croo, args.output, metafiles[0])
logger.info(cmd)
os.system(cmd)

bamfiles = glob.glob(args.output + "/**/*.bam", recursive = True)
for bamfile in bamfiles:
  cmd = "samtools index %s " % bamfile
  logger.info(cmd)
  os.system(cmd)

logger.info("done")