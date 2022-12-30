import os,sys,inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
cqsdir = os.path.abspath(os.path.dirname(currentdir) + "/CQS")
sys.path.insert(0,cqsdir) 

import argparse
import logging
import gzip
import subprocess
from collections import OrderedDict
from FileListUtils import readHashMap 
from Bio import SeqIO, bgzf

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Combine VCFs from GATK4/spin/trim",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Vcf files', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NotDEBUG)

args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/cqs/shengq2/jennifer_pietenpol_projects/20210502_brian_rnaseq_6151_hg38/gatk4_refine_SNV_04_merge/result/Breastspore_6151__mg_fileList1.list"
  args.output = "/scratch/cqs/shengq2/jennifer_pietenpol_projects/20210502_brian_rnaseq_6151_hg38/gatk4_refine_SNV_04_merge/result/Breastspore_6151.vcf.gz"
  #args.output = "/scratch/cqs/shengq2/jennifer_pietenpol_projects/20210502_brian_rnaseq_6151_hg38/gatk4_refine_SNV_04_merge/result/Breastspore_6151.vcf"

logger = logging.getLogger('combineVCFs')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

class VcfItem(object):
  def __init__(self, parts):
    self.CHROM = parts[0]
    self.POS = int(parts[1])
    self.ID = parts[2]
    self.REF = parts[3]
    self.ALT = parts[4]
    self.QUAL = float(parts[5])
    self.FILTER = parts[6]
    self.INFO = parts[7]
    self.FORMAT = parts[8]
    self.SNV = parts[9]

class VcfFile(object):
  def __init__(self):
    self.Headers = []
    self.Sample_names = []
    self.Format = ""
    self.Items = []
  
  def get_chroms(self):
    result = []
    for line in self.Headers:
      if line.startswith('##contig=<ID='):
        result.append(line.split('ID=', 1)[1].split(',',1)[0])
    #print(result)
    return(result)

  def read_pass_from_file(self, filepath):
    self.Items = []

    if filepath.endswith('.gz'):
      fh = gzip.open(filepath, "rt")
    else:
      fh = open(filepath, "rt") 
    with fh:
      for line in fh:
        if line.startswith('##'):
          self.Headers.append(line)
          continue

        if line.startswith('#'):
          self.Sample_names = line.rstrip().split('\t')[9:]
          break
      
      for line in fh:
        parts = line.rstrip().split('\t')
        if not parts[6] == "PASS":
          continue
        self.Items.append(VcfItem(parts))

def fill_valid_map(filename, filepath, valid_map):
  if filepath.endswith('.gz'):
    fh = gzip.open(filepath, "rt")
  else:
    fh = open(filepath, "rt") 
  with fh:
    for line in fh:
      if line.startswith('##'):
        continue
      if line.startswith('#'):
        break
    
    for line in fh:
      parts = line.rstrip().split('\t')
      
      chrom = parts[0]
      if chrom not in valid_map:
        continue
      cmap = valid_map[chrom]

      pos = int(parts[1])
      if not pos in cmap:
        continue
      amap = cmap[pos]
      
      allele_key =  parts[3] + "_" + parts[4]
      if allele_key in amap:
        values = amap[allele_key]
        if filename not in values[3]:
          values[3][filename] = parts[9]

fileMap = readHashMap(args.input)
#sample_names = [k for k in fileMap.keys()][0:2]
sample_names = [k for k in fileMap.keys()]

logger.info("reading valid SNVs ...")

vmap = {}
vcf1 = None
itemFormat = ""
for name in sample_names:
  filepath=fileMap[name][0]
  logger.info("  reading " + name + " ...")

  vf = VcfFile()
  vf.read_pass_from_file(filepath)

  if vf.Items.count == 0:
    continue

  for item in vf.Items:
    cmap = vmap.setdefault(item.CHROM, {})
    pmap = cmap.setdefault(item.POS,{})
    values = pmap.setdefault(item.REF + "_" + item.ALT, [item.REF, item.ALT, 0.0, {}])
    if item.QUAL > values[2]:
      values[2] = item.QUAL
    values[3][name] = item.SNV

  if vcf1 is None:
    vcf1 = vf
    itemFormat = vcf1.Items[0].FORMAT
    vf.Items.clear()

logger.info("filling valid SNVs ...")
for name in sample_names:
  filepath=fileMap[name][0]
  logger.info("  filling " + name + " ...")

  fill_valid_map(name, filepath, vmap)

chroms = vcf1.get_chroms()

with bgzf.BgzfWriter(args.output, "wb") as fout:
#with open(args.output, "wt") as fout:
  for line in vcf1.Headers:
    fout.write(line)

  fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % "\t".join(sample_names))
  for chrom in chroms:
    logger.info("Processing " + chrom + " ...")
    if not chrom in vmap:
      continue
    cmap = vmap[chrom]
    for pos in sorted(cmap.keys()):
      pmap = cmap[pos]
      for ra in sorted(pmap.keys()):
        values = pmap[ra]
        if not values[3]: #no PASS
          continue

        fout.write(f"{chrom}\t{pos}\t.\t{values[0]}\t{values[1]}\t{values[2]}\tPASS\t\t{itemFormat}")

        smap = values[3]
        for sample in sample_names:
          if sample in smap:
            fout.write("\t" + smap[sample])
          else:
            fout.write("\t./.")
        fout.write("\n")

logger.info("done.")
