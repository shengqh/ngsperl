import argparse
import sys
import logging
import os
import gzip
import subprocess

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Combine gCNV from GATK4 cohort pipeline",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input gCNV files', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file file", required=NotDEBUG)
parser.add_argument('-b', '--bedfile', action='store', nargs='?', help="Interval file in bed format", required=NotDEBUG)
parser.add_argument('-s', '--minimumScoreDifference', action='store', nargs='?', help="The minimum phred-scaled log posterior score difference between CNV event and normal event", default=30)
parser.add_argument('-f', '--minimumDuplicationFold', action='store', nargs='?', help="The minimum copy number fold change as duplication", default=2)
parser.add_argument('-p', '--percentage', action='store', default=0.9, type=float, nargs='?', help='Max sample percentage allowed')
parser.add_argument('--annovar_db', action='store', nargs='?', help='Annovar database folder')
parser.add_argument('--annovar_buildver', action='store', nargs='?', help='Annovar genome buildver')

args = parser.parse_args()

if DEBUG:
  args.input = "T:/Shared/Labs/Linton Lab/20180913_linton_exomeseq_2118_human_cutadapt/GATK4_CNV_Germline_CombineGCNV/result/linton_exomeseq_2118__fileList1.list"
  args.output = "T:/Shared/Labs/Linton Lab/20180913_linton_exomeseq_2118_human_cutadapt/GATK4_CNV_Germline_CombineGCNV/result/linton_exomeseq_2118.tsv"
  args.bedfile = "T:/Shared/Labs/Linton Lab/20180913_linton_exomeseq_2118_human_cutadapt/xgen-exome-research-panel-targetsae255a1532796e2eaa53ff00001c1b3c.nochr.bed"

logger = logging.getLogger('combineGCNV')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

fileMap = {}
with open(args.input) as fh:
  for line in fh:
    filepath, name = line.strip().split('\t', 1)
    fileMap[name] = filepath.strip()
    
samples = []
vcf1 = []
vcfMap = {}
bFirst = True
for name in fileMap.keys():
  filepath=fileMap[name]
  logger.info("reading " + name + " ...")
  samples.append(name)
  with gzip.open(filepath, "rt") as fh:
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
  
logger.info("reading " + args.bedfile + " ...")
annotationMap = {}
with open(args.bedfile, "r") as fin:
  for line in fin:
    parts = line.rstrip().split('\t')
    if not parts[0] in annotationMap:
      chrList = []
      annotationMap[parts[0]] = chrList
    else:
      chrList = annotationMap[parts[0]]
    if(len(parts) >= 4):
      chrList.append([int(parts[1]), int(parts[2]), parts[3]])
    else:
      chrList.append([int(parts[1]), int(parts[2]), "%s:%s-%s" % (parts[0], parts[1], parts[2])])
  
samples = sorted(samples)  
# samples = sorted(fileMap.keys())

combinedFile = args.output + ".combined.txt"
annovarInputFile = args.output + ".avinput"
logger.info("outputing to " + combinedFile + " ...")
with open(annovarInputFile, "w") as fav:
  with open(combinedFile, "w") as fout:
    fout.write("#Chrom\tStart\tEnd\tName\tGene\tcytoBand\t%s\n" % "\t".join(samples))
    intervalCount = len(vcf1)
    for idx in range(0, intervalCount):
      chrom = vcf1[idx][0]
      values = []
      for sample in samples:
        gt = vcfMap[sample][idx]
        if (gt.startswith("0:")):
          values.append("")
        else:
          parts = gt.split(":")
          cn = int(parts[1])
          expectCN = 2
          if chrom == 'X' or chrom == 'Y':
            expectCN = 1
          scores = parts[2].split(",")
          cnScore = int(scores[cn])
          expectScore = int(scores[expectCN])
          diffScore = expectScore - cnScore
          if cn > expectCN:
            if cn < expectCN * args.minimumDuplicationFold:
              values.append("")
              continue
            cnType = "DUP"
          else:
            cnType = "DEL"
          
          if expectScore - cnScore < args.minimumScoreDifference:
            values.append("")
          else:
            values.append("%s,%s,%s,%d" % (cnType, parts[0], cn, diffScore))
          
      if (all(gt == "" for gt in values)):
        continue

      chrom = vcf1[idx][0]
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
      
      cnvs = [v for v in values if v != ""]
      if len(cnvs) < len(values) * args.percentage:
        fout.write("%s\t%d\t%d\t%s\t\t\t%s\n" % (chrom, start, end, annotation, "\t".join(values)))
        fav.write("%s\t%d\t%d\t0\t0\n" % (chrom, start, end))

if args.annovar_db is None:
  if os.path.isfile(args.output):
    os.remove(args.output)
  os.rename(combinedFile, args.output)
else:
  logger.info("performing annovar ...")
  annovarOutput = annovarInputFile + ".annovar"
  subprocess.call(['table_annovar.pl', annovarInputFile, args.annovar_db, '-buildver', args.annovar_buildver, '-protocol', 'refGene,cytoBand','-operation', 'g,r', '--remove', '--outfile', annovarOutput])

  annovarOutputFile = annovarOutput + ".%s_multianno.txt" % args.annovar_buildver
  with open(args.output, "wt") as fout:
    with open(combinedFile, "rt") as fin:
      fout.write(fin.readline())
      with open(annovarOutputFile, "rt") as fann:
        fann.readline()
        for line in fin:
          lineAnno = fann.readline()
          annoParts = lineAnno.split('\t')
          parts = line.split('\t')
          parts[4] = annoParts[6]
          parts[5] = annoParts[10].rstrip()
          fout.write("\t".join(parts))

  os.remove(annovarOutputFile)
  os.remove(combinedFile)

logger.info("done.")
