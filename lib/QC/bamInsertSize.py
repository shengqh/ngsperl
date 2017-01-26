import pysam
import argparse
import sys
import logging
import os
from asyncore import read

parser = argparse.ArgumentParser(description="Build insert size distribution in SAM/BAM file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file (use "-" as stdin)', required=True)
parser.add_argument('--min-mapq', action='store', nargs='?', type=int, default=10, help="Minimum mapping quality of read")
parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output mismatch file name")

args = parser.parse_args()

logger = logging.getLogger('bamInsertSize')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

if args.output == "-":
  output = sys.stdout
else:
  tmpfile = args.output + ".tmp"
  output = open(tmpfile, 'w')
  
output.write("File\tInsertSize\tReadCount\n")
try:
  processed = 0
  inputfiles = args.input.split(',')
  for inputfile in inputfiles:
    if inputfile.lower().endswith(".bam"):
      openmode = "rb"
    else:
      openmode = "r"
    
    filename = "STDIN" if inputfile == "-" else os.path.splitext(os.path.basename(inputfile))[0]

    logger.info("processing %s ..." % filename)
    
    samfile = pysam.Samfile(inputfile, openmode)
    
    try:
      lengthMap = {}
      for read in samfile.fetch(until_eof=True):
        processed += 1
        if processed % 1000000 == 0:
          logger.info("processed %d" % processed)
          
        if read.is_unmapped or read.mapq < args.min_mapq or read.is_secondary or read.is_qcfail or read.is_duplicate or read.is_supplementary:
          continue;
        
        insertSize = read.template_length
        if(lengthMap.has_key(insertSize)):
          lengthMap[insertSize] = lengthMap[insertSize] + 1
        else:
          lengthMap[insertSize] = 1
      
      logger.info("total processed %d" % processed)
      for key in sorted(lengthMap.keys()):
        output.write("%s\t%d\t%d\n" % (filename, key, lengthMap[key]))
      
    finally:
      samfile.close()
      
  if args.output != "-":
    if os.path.isfile(args.output):
      os.remove(args.output)
    os.rename(tmpfile, args.output)
finally:
  output.close()            
