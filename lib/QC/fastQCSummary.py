import pysam
import argparse
import sys
import logging
import os
from asyncore import read

def initialize_logger(logfile, logname, isdebug):
  logger = logging.getLogger(logname)
  loglevel = logging.DEBUG if isdebug else logging.INFO
  logger.setLevel(loglevel)

  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')    
 
  # create console handler and set level to info
  handler = logging.StreamHandler()
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  # create error file handler and set level to error
  handler = logging.FileHandler(logfile, "w")
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  return(logger)

parser = argparse.ArgumentParser(description="Summerize the fastqc result.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input fastqc data file list', required=True)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output summary prefix", required=False)

args = parser.parse_args()

logger = initialize_logger(args.output + ".log", 'fastQCSummary', False)

logger.info("Output overall summary file.")
with open(args.output + ".summary.txt", "w") as fout:
  fout.write("Status\tCategory\tRawFile\tSample\n")
  with open(args.input, "r") as flistin:
    for line in flistin:
      parts = line.split('\t')
      datafile = parts[0]
      datafolder, datafilename = os.path.split(datafile)
      summaryfile = datafolder + "/summary.txt"
      with open(summaryfile, "r") as fsummary:
        for sline in fsummary:
          fout.write(sline.rstrip() + "\t" + parts[1].rstrip() + "\n")
    

  

