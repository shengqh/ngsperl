import argparse
import sys
import logging
import os
import errno
from asyncore import read

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

class MutectItem:
  def __init__(self, sampleName, line, normalIndex, tumorIndex):
    parts = line.rstrip().split("\t")

    self.SampleName = sampleName
    self.CHROM = parts[0]
    self.POS = parts[1]
    self.ID = parts[2]
    self.REF = parts[3]
    self.ALT = parts[4]
    self.QUAL = parts[5]
    self.FILTER = parts[6]
    self.INFO = parts[7]
    self.FORMAT = parts[8]
    self.NormalData = parts[normalIndex]
    self.TumorData = parts[tumorIndex]
    self.LocusKey = "_".join(self.CHROM, self.POS, self.REF, self.ALT)

    self.initDepth()

  def findDepth(self, parts, DP_index):
    result = int(parts[DP_index])
    return(result)

  def findMinorAllele(self, parts, AD_index):
    ad = parts[AD_index].split(',')
    result = int(ad[1])
    return(result)

  def initDepth(self):
    formatParts = self.FORMAT.split(':')
    normalParts = self.NormalData.split(':')
    tumorParts = self.TumorData.split(':')

    DP_index = formatParts.index("DP")
    self.NormalDepth = self.findDepth(normalParts, DP_index)
    self.TumorDepth = self.findDepth(tumorParts, DP_index)

    AD_index = formatParts.index("AD")
    self.TumorMinorAllele = self.findMinorAllele(tumorParts, AD_index)

class MutectResult:
  def clear(self):
    self.Comments = []
    self.ChromosomeItemMap = {}
    self.NormalSampleName = ""
    self.TumorSampleName = ""

  def __init__(self):
    self.clear()

  def findChromosomeFromComments(self):
    contigKey = "##contig=<ID="
    result = []
    for comment in self.Comments:
      if comment.startswith(contigKey):
        nextPart = comment[len(contigKey):]
        parts = nextPart.split(',', 1)
        result.append(parts[0])
    return(result)

  def findSampleName(self, line, sampleKey):
    sampleKeyEQ = sampleKey + "="
    if not sampleKeyEQ in line:
      raise Exception("The file is not mutect format, I cannot find %s in GATKCommandLine: %s" % (sampleKey, line))
    parts = line.split(sampleKeyEQ)[1]
    return(parts.split(" ")[0])

  def readFromFile(self, fileName, filePath):
    self.clear()
    with open(filePath, "rt") as fin:
      self.TumorSampleName = ""
      self.NormalSampleName = ""
      for line in fin:
        if line.startswith("##"):
          if line.startswith("##GATKCommandLine"):
            self.TumorSampleName = self.findSampleName(line, "tumor_sample_name")
            self.NormalSampleName = self.findSampleName(line, "normal_sample_name")
          else:
            self.Comments.append(line.rstrip())
        elif line.startswith("#CHROM"):
          if self.TumorSampleName == "":
            raise Exception("The file is not mutect format, I cannot find ##GATKCommandLine in %s" % args.input)
          parts = line.rstrip().split("\t")
          tumorIndex = parts.index(self.TumorSampleName)
          normalIndex = parts.index(self.NormalSampleName)
          logger.info("file=%s; tumor=%s; tumor_index=%d" % (os.path.basename(fileName), self.TumorSampleName, tumorIndex))
        else:
          item = MutectItem(fileName, line, normalIndex, tumorIndex)
          self.ChromosomeItemMap.setdefault(item.CHROM, []).append(item)

def filterMutect(logger, listFile, outputFile):  
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
    mutect.readFromFile(filePath)
    fileValueMap[fileName] = mutect
    if len(chroms) == 0:
      chroms = mutect.findChromosomeFromComments()
      comments = mutect.Comments

  with open(outputFile, "wt") as fout:
    for comment in Comments:
      fout.write("%s\n" % comment)
    fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % "\t".join(fileNames))

    for chrom in chroms:
      items = []
      for mutect in finalValueMap.values():
        items.extends(mutect.ChromosomeItemMap[chrom])
      
      posMap = {}
      for item in items:
        posMap.setdefault(item.POS, {}).setdefault(item.LosusKey, {})[item.SampleName] = item

      for pos in sorted(posMap.keys()):
        locusMap = posMap[pos]
        for locus in sorted(locusMap.keys()):
          sampleMap = locusMap[locus]
          item = [sampleMap.values()][0]
          fout.write("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (item.CHROM, item.POS, item.ID, item.REF, item.ALT, item.QUAL, item.FILTER, item.INFO, item.FORMAT))
          for sampleName in sampleNames:
            if sampleName in sampleMap:
              fout.write("\t%s" % sampleMap[sampleName].TumorData)
            else:
              fout.write("\t./.")
          fout.write("\n")
    
def main():
  DEBUG=True
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

  filterMutect(logger, args.input, args.output)

  logger.info("done.")

if __name__ == "__main__":
    main()
