import argparse
import sys
import logging
import os
import csv

class ReadItem:
  def __init__(self, sequence):
    self.Sequence = sequence
    self.TotalCount = 0
    self.SampleMap = {}

def getTotalCount(value):
  return -value.TotalCount

def build(logger, input, output, top):
  shortFileMap = {}
  samples = []
  with open(input, 'r') as sr:
    for line in sr:
      parts = line.rstrip().split('\t')
      samples.append(parts[1])
      shortFileMap[parts[1]] = parts[0]

  shortReadMap = {}
  for sample in samples:
    sampleFile = shortFileMap[sample]

    logger.info("  Reading " + sampleFile + " ...")
    with open(sampleFile, 'r') as fin:
      fin.readline()
      for line in fin:
        reads = line.rstrip().split('\t')
        count = int(reads[1])
        shortSeq = reads[2]

        read = shortReadMap.setdefault(shortSeq, ReadItem(shortSeq))
        read.SampleMap[sample] = count
        read.TotalCount += count

  allReads = shortReadMap.values()

  allReads = sorted(allReads, key=getTotalCount)
  if len(allReads) > top:
    allReads = allReads[0:top]

  with open(output, "wt") as fout:
    fout.write("Read\t%s\n" % "\t".join(samples) )
    for read in allReads:
      fout.write(read.Sequence)
      for sample in samples:
        if sample in read.SampleMap:
          fout.write("\t{}".format(read.SampleMap[sample]))
        else:
          fout.write("\t0")
      fout.write("\n")

  logger.info("Done.")

def main():
  parser = argparse.ArgumentParser(description="Build short read table",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG=False
  NOT_DEBUG = not DEBUG
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input short read list file', required=NOT_DEBUG)
  parser.add_argument('-t', '--top', action='store', nargs='?', type=int, default=200, help='Keep top X reads only (default 200)')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output table file", required=NOT_DEBUG)
  
  if NOT_DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/vickers_lab/projects/20200625_4893_2_RA_smRNA_mouse_v5_byTiger/host_genome/short_reads_table/result/RA_4893_2__fileList1.list"
    args.output = "/scratch/vickers_lab/projects/20200625_4893_2_RA_smRNA_mouse_v5_byTiger/host_genome/short_reads_table/result/RA_4893_2.count.txt"
  
  logger = logging.getLogger('shortReadTable')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  build(logger, args.input, args.output, args.top)
  
if __name__ == "__main__":
    main()
