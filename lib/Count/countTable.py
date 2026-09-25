import argparse
import sys
import logging
import os
import os.path
import logging
import errno
import shutil
from collections import OrderedDict

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def readFileMap(fileName):
  # no header, first column is name, second column is filepath
  check_file_exists(fileName)

  result = OrderedDict()
  with open(fileName) as fh:
    for line in fh:
      line = line.rstrip('\n')
      if line == '':
        continue
      name, filepath = line.split('\t', 1)
      filepath = filepath.strip()
      check_file_exists(filepath)
      result[name.strip()] = filepath
  return(result)

def readGeneCounts(count_file):
  # featureCounts output: lines starting with '#' are comments, header row starts with "Geneid";
  # first column is Geneid, seventh column is the count
  gene_counts = OrderedDict()
  with open(count_file, "rt") as fin:
    for line in fin:
      if line.startswith('#'):
        continue
      parts = line.rstrip('\n').split('\t')
      if parts[0] == 'Geneid':
        continue
      gene_counts[parts[0]] = parts[6]
  return(gene_counts)

OUTPUT_META_COLUMNS = ["length", "chr", "start", "end", "gene_biotype", "gene_name"]

def readGeneMap(fileName):
  # header present, first column is gene_id, remaining columns are reordered to OUTPUT_META_COLUMNS
  check_file_exists(fileName)

  with open(fileName) as fh:
    header = fh.readline().rstrip('\n').split('\t')
    col_indice = {col: idx for idx, col in enumerate(header[1:])}
    meta_columns = ["Feature_%s" % col for col in OUTPUT_META_COLUMNS]

    gene_map = OrderedDict()
    for line in fh:
      parts = line.rstrip('\n').split('\t')
      gene_id = parts[0]
      values = parts[1:]
      gene_map[gene_id] = [values[col_indice[col]] for col in OUTPUT_META_COLUMNS]
  return(meta_columns, gene_map)

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Get gene count table",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count file list (name<TAB>filepath, no header)', required=NotDEBUG)
parser.add_argument('-m', '--name_map', action='store', nargs='?', help='Gene meta file with header (gene_id<TAB>gene_name<TAB>...)', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/nobackup/shah_lab/shengq2/20260918_Qiagen_RNAseq_Kits/genetable/pbs/Qiagen_EV_hg38_tb.filelist"
  args.name_map = "/data/cqs/references/gencode/GRCh38.p13/gencode.v43.annotation.gtf.map"
  args.output = "/nobackup/shah_lab/shengq2/20260918_Qiagen_RNAseq_Kits/genetable/result/Qiagen_EV_hg38.test.count"

logger = logging.getLogger('countTable')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

count_files = readFileMap(args.input)
counts = OrderedDict()

total_files = len(count_files)
for idx, (name, count_file) in enumerate(count_files.items(), start=1):
  logger.info("Reading %d/%d : %s" % (idx, total_files, count_file))
  gene_counts = readGeneCounts(count_file)
  for geneid, count in gene_counts.items():
    counts.setdefault(geneid, {})[name] = count

meta_columns, gene_map = readGeneMap(args.name_map)

sample_names = list(count_files.keys())

def writeCountTable(output_file, biotype_filter=None):
  with open(output_file, "wt") as fout:
    fout.write("%s\n" % "\t".join(["Feature"] + meta_columns + sample_names))
    for feature in sorted(counts.keys()):
      values = [counts[feature].get(name, "0") for name in sample_names]
      if all(float(v) == 0 for v in values):
        continue
      meta = gene_map.get(feature, [""] * len(meta_columns))
      if biotype_filter is not None and meta[OUTPUT_META_COLUMNS.index("gene_biotype")] != biotype_filter:
        continue
      fout.write("%s\n" % "\t".join([feature] + meta + values))

writeCountTable(args.output)

base, ext = os.path.splitext(args.output)
proteincoding_output = "%s.proteincoding%s" % (base, ext)
writeCountTable(proteincoding_output, biotype_filter="protein_coding")

logger.info("done.")
