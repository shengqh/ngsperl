import sys
import os
import logging
import argparse
import re
import xml.etree.ElementTree as ET
from CountXmlUtils import readCountXmlQueryMap

DEBUG = False

if DEBUG:
  inputFile="/scratch/cqs/shengq2/vickers/20181120_smallRNA_269_933_2002_human_as_parclip/t2c/gsnap_smallRNA_t2c/result/Control_media_Rep1.T2C.tsv.xml"
  countXml="/scratch/cqs/shengq2/vickers/20181120_smallRNA_269_933_2002_human_as_parclip/t2c/gsnap_smallRNA_count/result/Control_media_Rep1/Control_media_Rep1.count.mapped.xml"
  outputFile="/scratch/cqs/shengq2/vickers/20181120_smallRNA_269_933_2002_human_as_parclip/t2c/gsnap_smallRNA_t2c/result/Control_media_Rep1.T2C.tsv.xml.txt"
else:
  parser = argparse.ArgumentParser(description="Extract T2C reads file.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input t2c xml file')
  parser.add_argument('-c', '--countXml', action='store', nargs='?', help='Input count xml file')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output T2C read file")

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  countXml = args.countXml
  outputFile = args.output

logger = logging.getLogger('T2CReads')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def readDict(filename):
  result = {}
  with open(filename, "r") as sr:
    for line in sr:
      parts = line.split('\t')
      result[parts[1].strip()] = parts[0]
  return(result)

def writeSequence(sw, beginChar, beginPos, offset, lastPos, sequence):
  res = beginChar
  if beginPos < offset:
    res = ' ' * (offset-beginPos) + res
  res = (res + sequence).ljust(lastPos + 1)
  sw.write(res)
  
with open(outputFile, "w") as sw:
  sw.write('Index\tSequence\tName\tT2C\tTotal\n')  

  featureMap = {}
  tree = ET.parse(inputFile)
  root = tree.getroot()
  features = root.find('subjectResult')
  for featureGroup in features.findall('subjectGroup'):
    feature = featureGroup.find('subject')
    featureName = feature.get("name")
    if featureName.startswith("miRNA:"):
      featureMap[featureName] = 1
  
  tree = ET.parse(countXml)
  root = tree.getroot()

  queryMap = {}
  queries = root.find('queries')
  for query in queries.findall('query'):
    query_name = query.get("name")
    for loc in query.findall('location'):
      seqname = loc.get("seqname")
      start = loc.get("start")
      end = loc.get("end")
      strand = loc.get("strand")
      cigar = loc.get("cigar")
      qkey = "%s_%s:%s-%s:%s" %(query_name, seqname, start, end, strand)
      queryMap[qkey] = cigar

  index = 0
  features = root.find('subjectResult')
  for featureGroup in features.findall('subjectGroup'):
    feature = featureGroup.find('subject')
    featureName = feature.get("name")
    
    if not featureName in featureMap:
      continue
    
    index = index + 1
    region = feature.find('region')
    sequence = region.get("sequence")
    strand = region.get("strand")
    
    beginPos = 0
    endPos = 0
    qnamelen = 0
    totalCount = 0
    t2cCount = 0
    for query in region.findall('query'):
      offset = int(query.get('offset'))
      seq_len = int(query.get('seq_len'))
      last = offset + seq_len
      beginPos = min([offset, beginPos])
      endPos = max([endPos, last])
      qname = query.get('qname')
      qnamelen = max(qnamelen, len(qname))
      query_count = int(query.get('query_count'))
      nnpm = int(query.get("nnpm"))
      if nnpm > 0:
        t2cCount = t2cCount + query_count
      totalCount = totalCount + query_count  
    sw.write("%d\t" % index)
    writeSequence(sw, '>', beginPos, 0, endPos, sequence)
    sw.write('\t%s\t%d\t%d\n' % (featureName, t2cCount, totalCount))  

    for query in region.findall('query'):
      qname = query.get("qname")
      loc = query.get("loc")
      qkey = "%s_%s" % (qname, loc)
      offset = int(query.get('offset'))
      cigar = queryMap[qkey]
      query_count = query.get('query_count')

      sw.write("%d\t" % index)
      writeSequence(sw, ' ', beginPos, offset, endPos, cigar)
      
      nnpm = int(query.get("nnpm"))
      if nnpm > 0:
        t2c = "Yes"
      else:
        t2c = "No"
      sw.write('\t%s\t%s\t%s\n' % (qname.ljust(qnamelen), t2c, query_count))  
        
logger.info("Result has been saved to %s" % outputFile)
