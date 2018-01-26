import csv
import argparse
import logging
import math
from Birdseed import BirdseedItem, readBirdseed

def isBaseValid(base):
  return base != " " and base != "-"

def buildBirdseedTable(listFile, outputFile):
  fileMap = {}
  with open(listFile, "r") as ir:
    header = ir.readline()
    for line in ir:
      parts = line.rstrip().split('\t')
      fileMap[parts[0]] = parts[1]
  
  filenames = sorted(fileMap.keys())
  bseeds = []
  print("writing %s ..." % outputFile)
  with open(outputFile, "w") as os:
    with open("%s.samples" % outputFile, "w") as osSample:
      bFirst=True
      fileIndex = 0
      for filename in filenames:
        fileIndex = fileIndex + 1
        print("reading %d/%d:%s ..." % (fileIndex, len(filenames), filename))
        file = fileMap[filename].split(',')[0]
        seeds = readBirdseed(file)
        
        if bFirst:
          os.write("Sample")
          for seed in seeds:
            os.write("\t%s" % seed.SNP)
          os.write("\n")
          bFirst = False
          bseeds = seeds
          bseedsIndecies = set()
          for idx in range(0, len(bseeds)):
            item = bseeds[idx]
            if isBaseValid(item.BaseA) and isBaseValid(item.BaseB):
              continue
            else:
              bseedsIndecies.add(idx)
        
        os.write(filename)
        osSample.write("%s\n" % filename)
        
        for seed in seeds:
          os.write("\t%s" % ("%0.2f" % math.log(seed.RatioBA, 2)))
        os.write("\n")
        
        if bseeds == seeds:
          continue
        
        keys = [idx for idx in bseedsIndecies]
        for idx in keys:
          item = seeds[idx]
          bitem = bseeds[idx]
          if isBaseValid(item.BaseA):
            bitem.BaseA = item.BaseA
          if isBaseValid(item.BaseB):
            bitem.BaseB = item.BaseB
          if isBaseValid(bitem.BaseA) and isBaseValid(bitem.BaseB):
            bseedsIndecies.remove(idx)
      
      os.write('BaseA')
      for item in bseeds:
        os.write("\t%s" % item.BaseA)
      os.write('\n')
      
      os.write('BaseB')
      for item in bseeds:
        os.write("\t%s" % item.BaseB)
      os.write('\n')
       
def main():
  parser = argparse.ArgumentParser(description="build birdseed table",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = True
  NOT_DEBUG = not DEBUG

  parser.add_argument('-i', '--inputFile', action='store', nargs='?', help='Input file list file', required=NOT_DEBUG)
  parser.add_argument('-o', '--outputFile', action='store', nargs='?', help='Output file', required=NOT_DEBUG)
  
  args = parser.parse_args()

  if DEBUG:
    args.inputFile = "/scratch/cqs/shengq2/guoyan/20180125_TCGA_eqtl_ratio/brca.filelist"
    args.outputFile = "/scratch/cqs/shengq2/guoyan/20180125_TCGA_eqtl_ratio/brca_intensity.tsv"

  print(args)
  
  logger = logging.getLogger('birdseed_table')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  buildBirdseedTable(args.inputFile, args.outputFile)

if __name__ == "__main__":
    main()
