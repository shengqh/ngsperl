import sys
import gzip
import os
import logging
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from CountXmlUtils import readCountXmlFeatures
from Feature import FeatureItem, FeatureGroup

DEBUG = 0

if DEBUG:
  #inputFile="/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_library/bowtie1_tRNA_pm_count/result/Liver_SRBIKO_11/Liver_SRBIKO_11.bam.count.mapped.xml"
  inputFile="/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_library/bowtie1_tRNA_pm_table/result/tRNA_pm_KCV_3018_77_78_79.count.xml"
  outputFile="/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/data_visualization/nonhost_library_tRNA_PositionVis_anticodon/tRNA_pm_KCV_3018_77_78_79.position.txt"
  fastaFile="/scratch/cqs/shengq2/references/smallrna/v3/GtRNAdb2/bowtie_index_1.1.2/GtRNAdb2.20161214.mature.fa"
  speciesMapFile = "/scratch/cqs/shengq2/references/smallrna/v3/GtRNAdb2/GtRNAdb2.20161214.map"
  species="bacteria"
else:
  parser = argparse.ArgumentParser(description="Generate smallRNA NTA read for Fastq file.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input xml file')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output coverage file")
  parser.add_argument('-f', '--fasta', action='store', nargs='?', help="Fasta file")
  parser.add_argument('-m', '--speciesMap', action='store', nargs='?', help="Species map file")
  parser.add_argument('-s', '--species', action='store', nargs='?', help="Required species name")
  parser.set_defaults(ccaa=False)

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  outputFile = args.output
  fastaFile=args.fasta
  speciesMapFile = args.speciesMap
  species = args.species

logger = logging.getLogger('nonhostLibraryCoverage')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def readFeatures(speciesMapFile, species, logger):
  result = {}
  logger.info("Reading species map from %s ..." % speciesMapFile)
  with open(speciesMapFile, "r") as f:
    header = f.readline()
    for line in f:
      parts = line.rstrip().split('\t')
      if parts[2] == species:
        result[parts[0]] = ""
  return result

def getAnticodon(tRNAname):
  m = re.search('tRNA-(.+?)-(\w+)', tRNAname)
  if m:
    return m.group(1) + m.group(2)
  
  m = re.search('trna\d+-(...)(...)', tRNAname)
  if m:
    return m.group(1) + m.group(2)
  
  return ("")

features = readFeatures(speciesMapFile, species, logger)
logger.info("%d features from %s" %(len(features), species))

logger.info("Reading feature-query in %s ..." % inputFile)
mappedFeatures = readCountXmlFeatures(inputFile)
logger.info("%d features mapped" % len(mappedFeatures))

for featureGroup in mappedFeatures:
  featureGroup.Features = [f for f in featureGroup.Features if f.Name in features]

mappedFeatures = [ mf for mf in mappedFeatures if len(mf.Features) > 0]
logger.info("%d features from %s mapped" % (len(mappedFeatures), species))

queries =set()
for fg in mappedFeatures:
  for query in fg.Queries:
    queries.add(query)
logger.info("%d reads from %d unique queries mapped" % (sum(q.Count for q in queries), len(queries)))

logger.info("Building anticodon map ...")
anticodonMap = {}
for featureGroup in mappedFeatures:
  for feature in featureGroup.Features:
    anticodon = getAnticodon(feature.Name)
    if anticodon == "":
      raise Exception("Wrong anticodon pattern: %s\n" % feature.Name)
    else:
      if anticodon not in anticodonMap:
        newGroup = FeatureGroup()
        anticodonMap[anticodon] = newGroup
      else:
        newGroup = anticodonMap[anticodon]
      newGroup.Features.append(feature)
      for query in featureGroup.Queries:
        newGroup.Queries.add(query)

logger.info("%d anticodon built" % len(anticodonMap))

logger.info("Reading sequences in %s ..." % fastaFile)
fasta_sequences = SeqIO.parse(open(fastaFile),'fasta')
for fasta in fasta_sequences:
  if fasta.id in features:
    features[fasta.id] = str(fasta.seq)

logger.info("Filling feature sequence ...")
for anticodon in anticodonMap:
  featureGroup = anticodonMap[anticodon]
  for feature in featureGroup.Features:
    feature.Sequence = features[feature.Name].upper()

def addCoverage(sample, querySequence, featureGroup):
  result = False
  for feature in featureGroup.Features:
    pos = feature.Sequence.find(querySequence)
    if pos != -1:
      result = True
      startpercentage = int(round(pos * 100.0 / len(feature.Sequence)))
      endpercentage = int(round((pos + len(querySequence) * 100.0 / len(feature.Sequence))))
      if not sample in featureGroup.Coverage:
        featureGroup.Coverage[sample] = [0] * 101
        featureGroup.CoverageCount[sample] = 0
        
      featureGroup.CoverageCount[sample] = featureGroup.CoverageCount[sample] + query.Count
      
      coverage = featureGroup.Coverage[sample]
      for idx in range(startpercentage, endpercentage):
        coverage[idx] = coverage[idx] + query.Count
        
      break
  return(result)
  
logger.info("Calculating coverage ...")

anticodons = sorted(anticodonMap.keys())
for anticodon in anticodons:
  featureGroup = anticodonMap[anticodon]
  for query in featureGroup.Queries:
    if addCoverage(query.Sample, query.Sequence, featureGroup):
      continue
    reverseCompSeq = str(Seq(query.Sequence).reverse_complement())
    if not addCoverage(query.Sample, reverseCompSeq, featureGroup):
      raise Exception("Cannot find %s or reverseCompSeq in %s" %(query.Sequence, featureGroup.Features[0].Sequence))

logger.info("Write result ...")
with open(outputFile, "w") as sw:
  sw.write("File\tFeature\tStrand\tTotalCount\tPositionCount\tPosition\tPercentage\n")
  for anticodon in anticodons:
    featureGroup = anticodonMap[anticodon]
    samples = sorted(featureGroup.Coverage.keys())
    for sample in samples:
      coverage = featureGroup.Coverage[sample]
      coverageCount = featureGroup.CoverageCount[sample]
      for idx in range(1,101):
        sw.write("%s\t%s\t-\t%d\t%d\t%d\t%.2f\n" %(sample, anticodon, coverageCount, coverage[idx], idx, (coverage[idx] * 1.0 / coverageCount)))
  
logger.info("Result has been saved to %s" % outputFile)
