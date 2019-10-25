import sys
import os
import logging
import argparse
import xml.etree.ElementTree as ET
from CountXmlUtils import readCountXmlFeatures
from os.path import basename

DEBUG = True
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="Extract Start Position From Xml",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count xml file', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output bam file", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="W:/SequenceProjects/3018-KCV-77_78_79/20191014_smallRNA_3018-KCV-77_78_79_mouse_v4_liverWT_byTiger/host_genome/bowtie1_genome_1mm_NTA_smallRNA_count/result/Liver_WT_12/Liver_WT_12.count.mapped.xml"
  args.output="W:/SequenceProjects/3018-KCV-77_78_79/20191014_smallRNA_3018-KCV-77_78_79_mouse_v4_liverWT_byTiger/host_genome/start_position/Liver_WT_12.startpos"
  
print(args)
  
logger = logging.getLogger('extractStartPosition')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def acceptFeatureName(featureName):
  return(featureName in { "miRNA:mmu-miR-122-5p", "tRNA:tRNA-Glu-CTC-1-1", "rRNA:rRNADB_28S_gi|120444900|ref|NR_003279.1|", "snoRNA:Gm22489:ENSMUSG00000080463.1","snRNA:Gm22317:ENSMUSG00000096206.1", "misc_RNA:Rn7s1:ENSMUSG00000099021.1"})

fileName = basename(args.input)

logger.info("Processing " + args.input)    
with open(args.output, "wt") as fout:
  fout.write("File\tFeature\tStrand\tPosition\tTotalCount\tPositionCount\tPercentage\tTotalUniqueCount\tPositionUniqueCount\tUniquePercentage\n")
  
  tree = ET.parse(args.input)
  root = tree.getroot()
  features = root.find('subjectResult')
  for featureGroup in features.findall('subjectGroup'):
    for feature in featureGroup.findall('subject'):
      featureName = feature.get("name")

      if not acceptFeatureName(featureName):
        continue

      for region in feature.findall('region'):
        strand = region.get("strand")
        size = int(region.get("size"))
        
        totalCount = 0
        totalUniqueCount = 0
        positionMap = {}
        uniquePositionMap = {}
        positionMap[0] = 0
        positionMap[size-1] = 0
        uniquePositionMap[0] = 0
        uniquePositionMap[size-1] = 0
        
        
        for query in region.findall('query'):
          offset = int(query.get('offset'))
          queryCount = int(query.get('query_count'))
          
          totalCount = totalCount + queryCount
          totalUniqueCount = totalUniqueCount + 1
          
          if not offset in positionMap:
            positionMap[offset] = queryCount
          else:
            positionMap[offset] = positionMap[offset] + queryCount
          
          if not offset in uniquePositionMap:
            uniquePositionMap[offset] = 1
          else:
            uniquePositionMap[offset] = uniquePositionMap[offset] + 1
        
        for offset in sorted(positionMap.keys()):  
          positionCount = positionMap[offset] if offset in positionMap else 0
          uniquePositionCount = uniquePositionMap[offset] if offset in uniquePositionMap else 0
          
          fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n" % 
                 (fileName, featureName, strand, offset, 
                  totalCount, positionCount, positionCount / totalCount,
                  totalUniqueCount, uniquePositionCount, uniquePositionCount / totalUniqueCount ))
        
        break #check the first region only
      
  