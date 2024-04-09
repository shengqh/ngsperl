import operator
from os.path import basename
import xml.etree.ElementTree as ET
from QueryItem import QueryItem
from Feature import FeatureItem, FeatureGroup

def readCountXmlQueryNames(fileName):
  result = set()
  tree = ET.parse(fileName)
  root = tree.getroot()
  queries = root.find('queries')
  for query in queries.findall('query'):
    query_name = query.get("name")
    result.add(query_name)
  return(result)

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

def defaultAcceptFeatureName(featureName):
  return(True)

def readCountXmlFeatures(fileName, acceptFunc = defaultAcceptFeatureName):
  queryMap = readCountXmlQueryMap(fileName)
  
  result = []
  tree = ET.parse(fileName)
  root = tree.getroot()
  features = root.find('subjectResult')
  for featureGroup in features.findall('subjectGroup'):
    fgroup = FeatureGroup()
    bAccept = False
    for feature in featureGroup.findall('subject'):
      featureName = feature.get("name")
      if acceptFunc(featureName):
        bAccept = True
        fi = FeatureItem(featureName, "")
        fgroup.Features.append(fi)
        for region in feature.findall('region'):
          for query in region.findall('query'):
            fgroup.Queries.add(queryMap[query.get('qname')])
            fi.EndPoints.append([int(query.get('offset')) + int(query.get('seq_len')), int(query.get('query_count'))])
    
    if bAccept:
      result.append(fgroup)
      for query in featureGroup.findall('query'):
        fgroup.Queries.add(queryMap[query.get('qname')])
        
  return(result)


def readCountXmlQueriesInFeatures(fileName, acceptFunc = defaultAcceptFeatureName):
  queryMap = readCountXmlQueryMap(fileName)
  
  result = set()
  tree = ET.parse(fileName)
  root = tree.getroot()
  features = root.find('subjectResult')
  for featureGroup in features.findall('subjectGroup'):
    fgroup = FeatureGroup()
    bAccept = False
    for feature in featureGroup.findall('subject'):
      featureName = feature.get("name")
      if acceptFunc(featureName):
        for region in feature.findall('region'):
          for query in region.findall('query'):
            result.add(queryMap[query.get('qname')])
            
  return(result)

def readCountXmlQueryLocationInFeatures(fileName):
  result = {}
  tree = ET.parse(fileName)
  root = tree.getroot()
  features = root.find('subjectResult')
  for featureGroup in features.findall('subjectGroup'):
    fgroup = FeatureGroup()
    for feature in featureGroup.findall('subject'):
      for region in feature.findall('region'):
        for query in region.findall('query'):
          queryName = query.get('qname')
          queryLocation = query.get('loc')
          if not queryName in result:
            result[queryName] ={queryLocation}
          else:
            result[queryName].add(queryLocation)
  return(result)
