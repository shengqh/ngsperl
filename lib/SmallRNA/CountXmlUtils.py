import operator
import xml.etree.ElementTree as ET
from QueryItem import QueryItem

def readCountXmlQueries(fileName, minCount):
  result = []
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
    if query_count >= minCount:
      result.append(QueryItem(query_sequence, query_name, query_count))
  result.sort(key=operator.attrgetter('QueryCount'), reverse=True)
  return(result)
  