import operator
from QueryItem import QueryItem

def readDupCountQueries(fileName, minCount, sorted=True):
  result = []
  with open(fileName, "r") as sr:
    sr.readline()
    for line in sr:
      parts = line.rstrip().split('\t')
      query_name = parts[0]
      query_count = int(parts[1])
      query_sequence = parts[2]
      if query_count >= minCount:
        result.append(QueryItem(query_name, query_sequence, query_count))
  if sorted:
    result.sort(key=operator.attrgetter('Count'), reverse=True)
  return(result)
  