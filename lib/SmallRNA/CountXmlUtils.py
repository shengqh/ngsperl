import operator
from os.path import basename
import xml.etree.ElementTree as ET
from QueryItem import QueryItem
from Feature import FeatureItem, FeatureGroup

def readCountXmlQueryMap(fileName):
  result = {}
  defaultSample = basename(fileName).split('.')[0]
  tree = ET.parse(fileName)
  root = tree.getroot()
  queries = root.find('queries')
  for query in queries.findall('query'):
    query_count = int(query.get("count"))
    if "sequence" in query.attrib:
      query_sequence=query.get("sequence")
    else:
      query_sequence=query.get("seq")
    query_name = query.get("name")
    qi = QueryItem(query_name, query_sequence, query_count)
    if "sample" in query.attrib:
      qi.Sample = query.get("sample")
    else:
      qi.Sample = defaultSample
    result[query_name] = qi
  return(result)
  
def readCountXmlQueries(fileName, minCount):
  resultMap = readCountXmlQueryMap(fileName)
  
  result = resultMap.values()
  if minCount > 1:
    result = [v for v in result if v.Count >= minCount]
    
  result.sort(key=operator.attrgetter('Count'), reverse=True)
    
  return(result)
  
def readCountXmlFeatures(fileName):
  queryMap = readCountXmlQueryMap(fileName);
  
  result = []
  tree = ET.parse(fileName)
  root = tree.getroot()
  features = root.find('subjectResult')
  for featureGroup in features.findall('subjectGroup'):
    fgroup = FeatureGroup()
    result.append(fgroup)
    for feature in featureGroup.findall('subject'):
      fgroup.Features.append(FeatureItem(feature.get("name"), ""))
    for query in featureGroup.findall('query'):
      fgroup.Queries.add(queryMap[query.get('qname')])
        
  return(result)
