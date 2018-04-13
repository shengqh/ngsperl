import collections

def readFileMap(fileName):
  result = collections.OrderedDict()
  
  with open(fileName, 'r') as f:
    for line in f:
      parts = line.rstrip().split('\t')
      if len(parts) > 1:
        result[parts[1]] = parts[0]

  return(result)
