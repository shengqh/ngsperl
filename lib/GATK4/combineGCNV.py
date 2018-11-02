import argparse
import sys
import logging
import os
import gzip

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Combine gCNV from GATK4 cohort pipeline",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input gCNV files', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file name", required=NotDEBUG)
parser.add_argument('-s', '--minimumScoreDifference', action='store', nargs='?', help="The minimum phred-scaled log posterior score difference between CNV event and normal event", default=30)

args = parser.parse_args()

if DEBUG:
  args.input = "T:/Shared/Labs/Linton Lab/20180913_linton_exomeseq_2118_human_cutadapt/GATK4_CNV_Germline_CombineGCNV/result/linton_exomeseq_2118__fileList1.list"
  args.output = "T:/Shared/Labs/Linton Lab/20180913_linton_exomeseq_2118_human_cutadapt/GATK4_CNV_Germline_CombineGCNV/result/linton_exomeseq_2118.tsv"

logger = logging.getLogger('combineGCNV')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

fileMap = {}
with open(args.input) as fh:
  for line in fh:
    file, name = line.strip().split('\t', 1)
    fileMap[name] = file.strip()
    
samples = []
vcf1 = []
vcfMap = {}
bFirst = True
for name, file in fileMap.iteritems():
  print("reading " + name + " ...")
  samples.append(name)
  with gzip.open(file, "rb") as fh:
    lines = []
    for line in fh:
      if line.startswith('#'):
        continue
      parts = line.rstrip().split('\t')
      if bFirst:
        vcf1.append(parts)
      lines.append(parts[9])
    vcfMap[name] = lines
  bFirst = False
  #if len(samples) == 2:
  #  break
  
annotationMap = {}
with open("T:/Shared/Labs/Linton Lab/20180913_linton_exomeseq_2118_human_cutadapt/xgen-exome-research-panel-targetsae255a1532796e2eaa53ff00001c1b3c.nochr.bed", "r") as fin:
  for line in fin:
    parts = line.split('\t')
    if not parts[0] in annotationMap:
      chrList = []
      annotationMap[parts[0]] = chrList
    else:
      chrList = annotationMap[parts[0]]
    chrList.append([int(parts[1]), int(parts[2]), parts[3]])
  
samples = sorted(samples)  
# samples = sorted(fileMap.keys())
with open(args.output, "w") as fout:
  fout.write("Locus\tName\t%s\n" % "\t".join(samples))
  intervalCount = len(vcf1)
  for idx in range(0, intervalCount):
    chr = vcf1[idx][0]
    values = []
    for sample in samples:
      gt = vcfMap[sample][idx]
      if (gt.startswith("0:")):
        values.append("")
      else:
        parts = gt.split(":")
        cn = int(parts[1])
        expectCN = 2
        if chr == 'X' or chr == 'Y':
          expectCN = 1
        scores = parts[2].split(",")
        cnScore = int(scores[cn])
        expectScore = int(scores[expectCN])
        diffScore = expectScore - cnScore
        if cn > expectCN:
          cnType = "DUP"
        else:
          cnType = "DEL"
        if expectScore - cnScore < args.minimumScoreDifference:
          values.append("")
        else:
          values.append("%s,%s,%s,%d" % (cnType, parts[0], cn, diffScore))
        
    if (all(gt == "" for gt in values)):
      continue

    chr = vcf1[idx][0]
    start = int(vcf1[idx][1])
    end = int(vcf1[idx][7][4:])
    
    annotation = "";
    annList = annotationMap[vcf1[idx][0]]
    for idx, ann in enumerate(annList):
      if ann[1] < start:
        continue
      if ann[0] > end:
        if idx == 0:
          annotation = ann[2]
        elif (ann[0] - end) < (start - annList[idx-1][1]):
          annotation = ann[2]
        else:
          annotation = annList[idx-1][2]
      else:
        annotation = ann[2]
      break 
    
    fout.write("%s:%d-%d\t%s\t%s\n" % (chr, start, end, annotation, "\t".join(values)))
