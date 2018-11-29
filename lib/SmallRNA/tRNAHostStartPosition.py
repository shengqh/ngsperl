import sys
import os
import logging
import argparse
import re
import xml.etree.ElementTree as ET

DEBUG = False

if DEBUG:
  inputFile="/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/data_visualization/host_genome_tRNA_start_position_vis/result/KCV_3018_77_78_79__fileList1.list"
  outputFile="/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/data_visualization/host_genome_tRNA_start_position_vis/result/KCV_3018_77_78_79.tsv"
else:
  parser = argparse.ArgumentParser(description="Generate tRNA start position file.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input xml file list')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output start position file")

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  outputFile = args.output

logger = logging.getLogger('tRNAHostStartPosition')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def readDict(filename):
  result = {}
  with open(filename, "r") as sr:
    for line in sr:
      parts = line.split('\t')
      result[parts[1].strip()] = parts[0]
  return(result)

xmlfiles = readDict(inputFile)

with open(outputFile, "w") as sw:
  sw.write("File\tFeature\tStrand\tTotalCount\tPositionCount\tPosition\tPercentage\n")

  sampleFeatureMap = {}  
  for sampleName in sorted(xmlfiles.keys()):
    logger.info("processing %s" % (sampleName))
    xmlfile = xmlfiles[sampleName]
    
    sampleFeatureMap[sampleName] = {}
  
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    features = root.find('subjectResult')
    for featureGroup in features.findall('subjectGroup'):
      feature = featureGroup.find('subject')
      featureName = feature.get("name")
      if featureName.startswith("tRNA:"):
        sampleFeatureMap[sampleName][featureName] = True
        region = feature.find('region')
        depth = {}
        totalCount = 0
        for query in region.findall('query'):
          offset = int(query.get('offset'))
          query_count = int(query.get('query_count'))
          totalCount = totalCount + query_count
          if offset in depth:
            depth[offset] = depth[offset] + query_count
          else:
            depth[offset] = query_count
        
        for idx in sorted(depth.keys()):
          sw.write("%s\t%s\t*\t%d\t%d\t%d\t%.3f\n" % (sampleName, featureName, totalCount,depth[idx], idx, depth[idx] * 1.0 / totalCount ))
  
  featureNames = set([fname for sf in sampleFeatureMap.values() for fname in sf.keys()])
  for sampleName in sorted(sampleFeatureMap.keys()):
    for featureName in featureNames:
      if featureName not in sampleFeatureMap[sampleName]:
        sw.write("%s\t%s\t*\t0\t0\t0\t0\n" % (sampleName, featureName ))
  
#dir_path = os.path.dirname(os.path.realpath(__file__))
#coverageR = dir_path + "/coverage.R"
#logger.info("Generate heatmap ...")
#cmd = "R --vanilla -f " + coverageR + " --args " + outputFile + " " + os.path.dirname(os.path.realpath(outputFile)) + "/"
#print(cmd + "\n")
#os.system(cmd)
        
logger.info("Result has been saved to %s" % outputFile)
