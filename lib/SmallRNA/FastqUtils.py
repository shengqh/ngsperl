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
          resultMap[sequence].QueryCount = resultMap[sequence].QueryCount + 1
        else:
          resultMap[sequence] = QueryItem(sequence, header[1:], 1)
  finally:
    f.close()
  
  result = resultMap.values()
  if minCount > 1:
    result = [v for v in result if v.QueryCount >= minCount]
    
  result.sort(key=operator.attrgetter('QueryCount'), reverse=True)
    
  return(result)
  