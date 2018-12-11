import sys
import os
import logging
import argparse
import re
import xml.etree.ElementTree as ET



def readDict(filename):
  result = {}
  with open(filename, "r") as sr:
    for line in sr:
      parts = line.split('\t')
      result[parts[1].strip()] = parts[0]
  return(result)

def extract(logger, inputFile, outputFile):
  xmlfiles = readDict(inputFile)
  
  with open(outputFile, "w") as sw:
    sw.write("File\tFeature\tStrand\tTotalCount\tPositionCount\tPosition\tPercentage\n")
  
    sampleFeatureMap = {}  
    for sampleName in sorted(xmlfiles.keys()):
      logger.info("processing %s" % (sampleName))
      xmlfile = xmlfiles[sampleName]
      
      sampleFeatureMap[sampleName] = {}
      
      with open(xmlfile, "r") as fin:
        for line in fin:
          if "<subject name=\"tRNA:" in line:
            xmlstring = line
            insubject = True
            for subline in fin:
              xmlstring = xmlstring + subline
              if "</subject" in subline:
                break
            
            feature = ET.fromstring(xmlstring)
            featureName = feature.get("name")
            sampleFeatureMap[sampleName][featureName] = True
            depth = {}
            totalCount = 0
            queries = {}
            for region in feature.findall('region'):
              for query in region.findall('query'):
                query_name = query.get('qname')
                if query_name in queries:
                  continue
                queries[query_name] = True
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


def main():
  DEBUG = True
  NOT_DEBUG = not DEBUG

  parser = argparse.ArgumentParser(description="Generate tRNA start position file.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input xml file list', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output start position file", required=NOT_DEBUG)

  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  print(args)

  if DEBUG:
    args.input="E:/sqh/programs/csharp/OmicsLabCSharp/CQS.Test/data/trna_start_position_fileList1.list"
    args.output="H:/shengquanhu/projects/temp/trna_start_position.tsv"

  logger = logging.getLogger('tRNAHostStartPosition')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  extract(logger, args.input, args.output)
  
if __name__ == "__main__":
    main()
