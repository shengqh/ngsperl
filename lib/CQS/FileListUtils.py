import sys
import os
import string
from collections import OrderedDict

def readHashMap(fileName):
  """Read file list file and return key/array of value hash map. Can be used to read group definition.
  
  Arguments:
    fileName {str} -- file name, in which the second column is key and first column is value, no header. One key may correspond to multiple values, such as group:file definition.
  
  Result:
    Hash map, key:array of values
  """
  result = OrderedDict()
  with open(fileName, "r") as f:
    for line in f:
      parts = line.rstrip().split('\t')
      if(len(parts) > 1):
        if parts[1] not in result:
          result[parts[1]] = [parts[0]]
        else:
          result[parts[1]].append(parts[0])
  return(result)

def readUniqueHashMap(fileName, keyIndex = 1):
  """Read file list file and return key/value hash map. Can be used to sample name/sample file definition, or parameter options definition.
  
  Arguments:
    fileName {str} -- the file containing at least two columns, one is key and another is value, no header. One key corresponds to one value.
    keyIndex {int} -- the column index (0-based) of key.
  
  Result:
    Hash map, key:value
  """

  valueIndex = 1 if keyIndex == 0 else 0
  maxIndex = max(valueIndex, keyIndex)

  result = OrderedDict()
  with open(fileName, "r") as f:
    for line in f:
      parts = line.rstrip().split('\t')
      if(len(parts) > maxIndex):
        result[parts[keyIndex]] = parts[valueIndex]

  return(result)

