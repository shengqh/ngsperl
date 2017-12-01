import operator
import gzip
from QueryItem import QueryItem

def readFastqQueries(fileName, minCount):
  resultMap = {}
  
  gzipped = fileName.endswith(".gz")
  if gzipped:
    f = gzip.open(fileName, 'rt')
  else:
    f = open(fileName, 'r')
  
  try:
    while True:
      header = f.readline()
      if '' == header:
        break

      if not header.startswith("@"):
        continue

      sequence = f.readline().strip()
      f.readline()
      f.readline()

      if len(sequence) > 0:
        if sequence in resultMap:
          resultMap[sequence].Count = resultMap[sequence].Count + 1
        else:
          resultMap[sequence] = QueryItem(header[1:], sequence, 1)
  finally:
    f.close()

  result = resultMap.values()
  if minCount > 1:
    result = [v for v in result if v.Count >= minCount]
    
  result.sort(key=operator.attrgetter('Count'), reverse=True)
    
  return(result)
  
def readFastqQueryMap(fileName):
  result = {}

  resultList = readFastqQueries(fileName, 1)
  for query in resultList:
    result[query.Name] = query;
  
  return (result)
