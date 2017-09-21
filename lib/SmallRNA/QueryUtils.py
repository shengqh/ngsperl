from CountXmlUtils import readCountXmlQueries
from DupCountUtils import readDupCountQueries
from FastqUtils import readFastqQueries

def readQueries(fileName, minCount):
  result = {}
  if fileName.endswith(".xml"):
    return readCountXmlQueries(fileName, minCount)
  elif fileName.endswith(".dupcount"):
    return readDupCountQueries(fileName, minCount)
  elif fileName.endswith(".fastq.gz") or fileName.endswith(".fastq.short.gz"):
    return readFastqQueries(fileName, minCount)
  else:
    raise Exception("Unknow file format %s" % fileName)
  