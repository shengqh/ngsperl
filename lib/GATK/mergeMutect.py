import argparse
import sys
import logging
import os
import errno
import gzip
from asyncore import read
from Mutect import MutectItem, MutectResult

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

def mergeMutect(logger, listFile, outputFile):  
  fileMap = readFileMap(listFile)
  checkFileMap(fileMap)

  fileValueMap = {}
  chroms = []
  comments = []

  fileNames = sorted(fileMap.keys())
  for fileName in fileNames:
    filePath = fileMap[fileName]

    logger.info("Reading %s ..." % filePath)

    mutect = MutectResult()
    mutect.readFromFile(logger, fileName, filePath)
    fileValueMap[fileName] = mutect
    if len(chroms) == 0:
      chroms = mutect.findChromosomeFromComments()
      comments = mutect.Comments

  has_normal = any(v.NormalSampleName != None for v in fileValueMap.values())

  logger.info("Output result to %s ..." % outputFile)
  with open(outputFile, "wt") as fout:
    for comment in comments:
      if comment.startswith("##INFO=<ID=LOD"):
        if has_normal:
          fout.write('##FORMAT=<ID=ND,Number=1,Type=Integer,Description="Approximate normal sample read depth (reads with MQ=255 or with bad mates are filtered)">\n')
        fout.write("%s\n" % comment.replace("##INFO=", "##FORMAT="))
      else:
        fout.write("%s\n" % comment)
    fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % "\t".join(fileNames))

    for chrom in chroms:
      items = []
      for mutect in fileValueMap.values():
        if chrom in mutect.ChromosomeItemMap:
          items.extend(mutect.ChromosomeItemMap[chrom])
      
      posMap = {}
      for item in items:
        posMap.setdefault(item.POS, {}).setdefault(item.LocusKey, {})[item.SampleName] = item

      for pos in sorted(posMap.keys()):
        locusMap = posMap[pos]
        for locus in sorted(locusMap.keys()):
          sampleMap = locusMap[locus]
          item = [v for v in sampleMap.values()][0]
          if has_normal:
            fout.write("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s:ND:LOD" % (item.CHROM, item.POS, item.ID, item.REF, item.ALT, item.QUAL, item.FILTER, item.INFO, item.FORMAT))
          else:
            fout.write("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s:LOD" % (item.CHROM, item.POS, item.ID, item.REF, item.ALT, item.QUAL, item.FILTER, item.INFO, item.FORMAT))

          for sampleName in fileNames:
            if sampleName in sampleMap:
              item = sampleMap[sampleName]
              if has_normal:
                fout.write("\t%s:%d:%s" % (item.TumorData, item.NormalDepth, item.LOD))
              else:
                fout.write("\t%s:%s" % (item.TumorData, item.LOD))
            else:
              fout.write("\t./.")
          fout.write("\n")
    
def main():
  DEBUG=False
  NotDEBUG=not DEBUG

  parser = argparse.ArgumentParser(description="merge mutect result and keep tumor sample only.",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input vcf list file', required=NotDEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output vcf file", required=NotDEBUG)

  args = parser.parse_args()

  if DEBUG:
    args.input = "H:/shengquanhu/projects/20190610_Ciombior_ExomeSeq/Ciombor_ExomeSeq__fileList1.list"
    args.output = "H:/shengquanhu/projects/20190610_Ciombior_ExomeSeq/combined.tumor.vcf"

  logger = logging.getLogger('mergeMutect')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  mergeMutect(logger, args.input, args.output)

  logger.info("done.")

if __name__ == "__main__":
    main()
