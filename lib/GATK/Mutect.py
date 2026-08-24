import gzip
import os
from collections import OrderedDict
from scipy.stats import fisher_exact

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
    self.Line = line

    if normalIndex != -1:
      self.NormalData = parts[normalIndex]
    else:
      self.NormalData = None

    self.TumorData = parts[tumorIndex]
    self.LocusKey = "_".join([self.CHROM, str(self.POS), self.REF, self.ALT])

    self.parseData()

  def findDepth(self, parts, DP_index):
    result = int(parts[DP_index])
    return(result)

  def findAlleleDepth(self, parts, AD_index):
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##The first value is the depth of the reference allele, and the second value is the depth of the alternate allele. 
    ad = parts[AD_index].split(',')
    result = [int(ad[0]), int(ad[1])]
    return(result)

  def findAltAlleleFrequency(self, parts, AF_index):
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">
    af = parts[AF_index].split(',')
    result = float(af[0])
    return(result)

  def parseData(self):
    formatParts = self.FORMAT.split(':')
    tumorParts = self.TumorData.split(':')

    DP_index = formatParts.index("DP")
    self.TumorDepth = self.findDepth(tumorParts, DP_index)

    if self.NormalData != None:
      normalParts = self.NormalData.split(':')
      self.NormalDepth = self.findDepth(normalParts, DP_index)

    AD_index = formatParts.index("AD")
    self.TumorRefAlleleDepth, self.TumorAltAlleleDepth = self.findAlleleDepth(tumorParts, AD_index)

    AF_index = formatParts.index("AF")
    self.TumorAltAlleleFrequency = self.findAltAlleleFrequency(tumorParts, AF_index)
    if self.NormalData != None:
      self.NormalRefAlleleDepth, self.NormalAltAlleleDepth = self.findAlleleDepth(normalParts, AD_index)
      self.NormalAltAlleleFrequency = self.findAltAlleleFrequency(normalParts, AF_index)
    
    self.FORMAT = self.FORMAT
    if ";LOD=" in self.INFO:
      lodParts = self.INFO.split(";LOD=")
      self.LOD = lodParts[1]
      self.INFO = lodParts[0]
    elif ";TLOD=" in self.INFO:
      lodParts = self.INFO.split(";TLOD=")
      self.LOD = lodParts[1]
      self.INFO = lodParts[0]
    else:
      self.LOD = None

  def getFisherPValue(self):
    if self.NormalData is None:
      return None

    ref_alt_table = [[self.NormalRefAlleleDepth, self.NormalAltAlleleDepth],
                     [self.TumorRefAlleleDepth, self.TumorAltAlleleDepth]]

    _, p_value = fisher_exact(ref_alt_table)
    return p_value

class MutectResult:
  def clear(self):
    self.Comments = []
    self.ChromosomeItemMap = OrderedDict()
    self.NormalSampleName = None
    self.TumorSampleName = ""
    self.Line = None

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

  def doFindSampleName(self, line, sampleKey, eqValue, name, required):
    sampleKeyEQ = sampleKey + eqValue
    if not sampleKeyEQ in line:
      if required:
        raise Exception("The file is not %s format, I cannot find %s in GATKCommandLine: %s" % (name, sampleKey, line))
      else:
        return(None)
    parts = line.split(sampleKeyEQ)[1]
    return(parts.split(" ")[0])

  def findMutect1SampleName(self, line, sampleKey, required):
    return self.doFindSampleName(line, sampleKey, "=", "mutect1", required)

  def findMutect2SampleName(self, line, sampleKey, required):
    return self.doFindSampleName(line, sampleKey, " ", "mutect2", required)

  def readFromFile(self, logger, fileName, filePath):
    self.clear()

    if filePath.endswith(".gz"):
      fin = gzip.open(filePath,'rt')
    else:
      fin = open(filePath, "rt")

    normalIndex = -1

    with fin:
      self.TumorSampleName = None
      self.NormalSampleName = None
      for line in fin:
        if line.startswith("##"):
          if line.startswith("##GATKCommandLine") and (self.TumorSampleName is None):
            #print(line)
            if "ID=Mutect2" in line:
              self.TumorSampleName = self.findMutect2SampleName(line, "--tumor-sample", True)
              self.NormalSampleName = self.findMutect2SampleName(line, "--normal-sample", False)
            elif "ID=MuTect" in line:
              self.TumorSampleName = self.findMutect1SampleName(line, "tumor_sample_name", True)
              self.NormalSampleName = self.findMutect1SampleName(line, "normal_sample_name", False)
          if line.startswith("##normal_sample="):
            self.NormalSampleName = line.rstrip().split('=')[1]
          if line.startswith("##tumor_sample="):
            self.TumorSampleName = line.rstrip().split('=')[1]
          self.Comments.append(line.rstrip())
        elif line.startswith("#CHROM"):
          if self.TumorSampleName == "":
            raise Exception("The file is not mutect format, I cannot find ##GATKCommandLine in %s" % filePath)
          parts = line.rstrip().split("\t")
          tumorIndex = parts.index(self.TumorSampleName)
          if self.NormalSampleName != None:
            normalIndex = parts.index(self.NormalSampleName)
          logger.info("file=%s; tumor=%s; tumor_index=%d" % (os.path.basename(fileName), self.TumorSampleName, tumorIndex))
        else:
          item = MutectItem(fileName, line, normalIndex, tumorIndex)
          if item.FILTER == "PASS":
            self.ChromosomeItemMap.setdefault(item.CHROM, []).append(item)

class MutectSampleItem:
  def __init__(self, data, normalDepth, tumorDepth, refAlleleDepth, altAlleleDepth):
    self.Data = data
    self.NormalDepth = normalDepth
    self.TumorDepth = tumorDepth
    self.RefAlleleDepth = refAlleleDepth
    self.AltAlleleDepth = altAlleleDepth

class MultiMutectItem:
  def __init__(self, line):
    self.Line = line
    self.Samples = []
    parts = line.rstrip().split('\t')
    formatParts = parts[8].split(':')
    if ":ND" in parts[8]:
      if formatParts[-2] == "ND":
        ND_index = -2
      else:
        ND_index = formatParts.index("ND")
    else:
      ND_index = None
    DP_index = formatParts.index("DP")
    AD_index = formatParts.index("AD")

    for sampleIndex in range(9, len(parts)):
      sampleData = parts[sampleIndex]
      if sampleData.startswith("./."):
        self.Samples.append(MutectSampleItem(sampleData, 0, 0, 0, 0))
        continue

      sampleParts = sampleData.split(':')

      if ND_index != None:
        normalDepth = int(sampleParts[ND_index])
      else:
        normalDepth = None
      tumorDepth = int(sampleParts[DP_index])
      alleleParts = sampleParts[AD_index].split(',')
      refAlleleDepth = int(alleleParts[0])
      altAlleleDepth = int(alleleParts[1])

      self.Samples.append(MutectSampleItem(sampleData, normalDepth, tumorDepth, refAlleleDepth, altAlleleDepth))
