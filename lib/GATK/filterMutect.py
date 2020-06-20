import argparse
import sys
import logging
import os
import errno
from asyncore import read

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="filter mutect result to keep tumor sample only.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input vcf list file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output vcf file", required=NotDEBUG)

args = parser.parse_args()

if DEBUG:
  args.input = "H:/shengquanhu/projects/20190610_Ciombior_ExomeSeq/Ciombor_ExomeSeq__fileList1.list"
  args.output = "H:/shengquanhu/projects/20190610_Ciombior_ExomeSeq/combined.tumor.vcf"

logger = logging.getLogger('filterMutect')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("Reading file list %s ..." % args.input)

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def readFileMap(fileName):
  check_file_exists(fileName)

  result = {}
  with open(fileName) as fh:
    for line in fh:
      filepath, name = line.strip().split('\t', 1)
      result[name] = filepath.strip()
  return(result)

def checkFileMap(fileMap):
  for sname in fileMap.keys():
    sfile = fileMap[sname]
    check_file_exists(sfile)

def readVcfFile(fileName):
  comments = []
  items = []
  header = ""
  with open(fileName, "rt") as fin:
    tumor_sample_name = ""
    for line in fin:
      if line.startswith("##"):
        if line.startswith("##GATKCommandLine"):
          if not "tumor_sample_name=" in line:
            raise Exception("The file is not mutect format, I cannot find tumor_sample_name in GATKCommandLine: %s" % line);
          tumor_parts = line.split("tumor_sample_name=")[1]
          tumor_sample_name = tumor_parts.split(" ")[0]
        else:
          comments.append(line.rstrip())
      elif line.startswith("#CHROM"):
        if tumor_sample_name == "":
          raise Exception("The file is not mutect format, I cannot find ##GATKCommandLine in %s" % args.input);
        parts = line.rstrip().split("\t")
        format_index = parts.index("FORMAT")
        tumor_index = parts.index(tumor_sample_name)
        logger.info("format_index=%d; tumor_index=%d" % (format_index, tumor_index))
        indecies = [index for index in range(0, len(parts)) if index <= format_index]
        header = "\t".join([parts[index] for index in indecies])
      else:
        parts = line.rstrip().split("\t")
        locus = [parts[index] for index in indecies]
        locus.append(parts[tumor_index])
        items.append(locus)
  return(items, comments, header)
  
fileMap = readFileMap(args.input)
checkFileMap(fileMap)

fileValueMap = {}
comments = []
header = ""

locusKeyMap = {}
fileNames = sorted(fileMap.keys())
for fileName in fileNames:
  filePath = fileMap[fileName]
  items, comments, header = readVcfFile(filePath)
  itemMap = {}
  for item in items:
    locusKey = "_".join([item[0], item[1], item[3], item[4]])
    itemMap[locusKey] = item
    if locusKey not in locusKeyMap:
      locusKeyMap[locusKey] = item
  fileValueMap[fileName] = itemMap

with open(args.output, "wt") as fout:
  for comment in comments:
    fout.write("%s\n" % comment)
  fout.write("%s\t%s\n" % (header, "\t".join(fileNames)))
  for locusKey in locusKeyMap:
    item = locusKeyMap[locusKey]
    fout.write("%s\tSOMATIC;VT=SNP\t%s" % ("\t".join(item[0:6]), item[8]))
    for fileName in fileNames:
      itemMap = fileValueMap[fileName]
      if locusKey in itemMap:
        locusItem = itemMap[locusKey]
        fout.write("\t%s" % locusItem[-1])
      else:
        fout.write("\t./.")
    fout.write("\n")
          
logger.info("done.")     
  