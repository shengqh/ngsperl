import gzip
import os

class MutectItem:
  def __init__(self, sampleName, line, normalIndex, tumorIndex):
    parts = line.rstrip().split("\t")

    self.SampleName = sampleName
    self.CHROM = parts[0]
    self.POS = int(parts[1])
    self.ID = parts[2]
    self.REF = parts[3]
    self.ALT = parts[4]
    self.QUAL = parts[5]
    self.FILTER = parts[6]
    self.INFO = parts[7]
    self.FORMAT = parts[8]
    self.LOD = 0
    self.NormalData = parts[normalIndex]
    self.TumorData = parts[tumorIndex]
    self.LocusKey = "_".join([self.CHROM, str(self.POS), self.REF, self.ALT])

    self.parseData()

  def findDepth(self, parts, DP_index):
    result = int(parts[DP_index])
    return(result)

  def findMinorAllele(self, parts, AD_index):
    ad = parts[AD_index].split(',')
    result = int(ad[1])
    return(result)

  def parseData(self):
    formatParts = self.FORMAT.split(':')
    normalParts = self.NormalData.split(':')
    tumorParts = self.TumorData.split(':')

    DP_index = formatParts.index("DP")
    self.NormalDepth = self.findDepth(normalParts, DP_index)
    self.TumorDepth = self.findDepth(tumorParts, DP_index)

    AD_index = formatParts.index("AD")
    self.TumorMinorAllele = self.findMinorAllele(tumorParts, AD_index)
    
    self.FORMAT = self.FORMAT
    lodParts = self.INFO.split(";LOD=")
    self.LOD = lodParts[1]
    self.INFO = lodParts[0]

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

  def readFromFile(self, logger, fileName, filePath):
    self.clear()

    if filePath.endswith(".gz"):
      fin = gzip.open(filePath,'rt')
    else:
      fin = open(filePath, "rt")

    with fin:
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

class MutectSampleItem:
  def __init__(self, data, normalDepth, tumorDepth, majorAlleleDepth, minorAlleleDepth):
    self.Data = data
    self.NormalDepth = normalDepth
    self.TumorDepth = tumorDepth
    self.MajorAlleleDepth = majorAlleleDepth
    self.MinorAlleleDepth = minorAlleleDepth

class MultiMutectItem:
  def __init__(self, line):
    self.Line = line
    self.Samples = []
    parts = line.rstrip().split('\t')
    formatParts = parts[8].split(':')
    ND_index = formatParts.index("ND")
    DP_index = formatParts.index("DP")
    AD_index = formatParts.index("AD")

    for sampleIndex in range(9, len(parts)):
      sampleData = parts[sampleIndex]
      if sampleData.startswith("./."):
        self.Samples.append(MutectSampleItem(sampleData, 0, 0, 0, 0))
        continue

      sampleParts = sampleData.split(':')
      normalDepth = int(sampleParts[ND_index])
      tumorDepth = int(sampleParts[DP_index])
      alleleParts = sampleParts[AD_index].split(',')
      majorAlleleDepth = int(alleleParts[0])
      minorAlleleDepth = int(alleleParts[1])

      self.Samples.append(MutectSampleItem(sampleData, normalDepth, tumorDepth, majorAlleleDepth, minorAlleleDepth))

class MultiMutectResult:
  def clear(self):
    self.Comments = []
    self.Data = []

  def __init__(self):
    self.clear()

  def readFromFile(self, filePath):
    self.clear()

    if filePath.endswith(".gz"):
      fin = gzip.open(filePath,'rt')
    else:
      fin = open(filePath, "rt")

    with fin:
      for line in fin:
        if line.startswith("#"):
          self.Comments.append(line)
        else:
          self.Data.append(MultiMutectItem(line))
