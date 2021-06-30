import argparse
import sys
import logging
import os

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

parser = argparse.ArgumentParser(description="Summarize the fastqc result.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input fastqc data file list', required=True)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output summary prefix", required=False)

args = parser.parse_args()

logger = initialize_logger(args.output + ".log", 'fastQCSummary', False)

logger.info("Output overall summary file.")

class Section:
  def __init__(self, fileName, sectionTag, headerTag):
    self.fileName = fileName
    self.sectionTag = sectionTag
    self.headerTag = headerTag
    self.inSection = False
    self.header = ""
    self.data = {}

sections = [Section("adapter", ">>Adapter Content", "#Position"), 
            Section("baseQuality", ">>Per base sequence quality", "#Base"),
            Section("baseContent", ">>Per base sequence content", "#Base"),
            Section("sequenceGC", ">>Per sequence GC content", "#GC Content")]
summary={}
reads={}
overrepresented={}
overrepresentedHeader = ""
has_error = False
with open(args.input, "r") as flistin:
  for line in flistin:
    parts = line.split('\t')
    sampleName = parts[1].rstrip()
    
    datafile = parts[0]
    datafolder, datafilename = os.path.split(datafile)
    datafilePrefix = os.path.splitext(os.path.basename(datafolder))[0]
    if datafilePrefix.endswith("_fastqc"):
      datafilePrefix = datafilePrefix[0:-7]
    
    if not sampleName in summary:
      summary[sampleName] = {}
    prefixDic = {}
    summary[sampleName][datafilePrefix] = prefixDic

    summaryfile = datafolder + "/summary.txt"
    if (not os.path.isfile(summaryfile)):
      logger.error("File not exists: %s" % summaryfile)
      has_error = True
      continue

    with open(summaryfile, "r") as fsummary:
      for sline in fsummary:
        sparts = sline.rstrip().split('\t')
        prefixDic[sparts[1]] = sparts[0]

    if not sampleName in reads:
      reads[sampleName] = {}

    if not sampleName in overrepresented:
      overrepresented[sampleName] = {}

    for section in sections:
      section.inSection = False
      section.done = False
      if not sampleName in section.data:
        section.data[sampleName] = {}      

    with open(datafile, "r") as fdata:
      bInOver = False
      bInAdapter = False
      for sline in fdata:
        if sline.startswith("Total Sequences"):
          reads[sampleName][datafilePrefix] = sline[16:].rstrip()
          continue

        if sline.startswith(">>Overrepresented"):
          bInOver = True
          continue

        if bInOver:
          if sline.startswith("#Sequence"):
            overrepresentedHeader = sline[1:].rstrip()
            continue
          if not sline.startswith(">>END_MODULE"):
            overrepresented[sampleName][datafilePrefix] = sline.rstrip()
          bInOver = False
          continue

        for section in sections:
          if sline.startswith(section.sectionTag):
            section.inSection = True
            section.data[sampleName][datafilePrefix] = []
            continue
          
          if section.inSection:
            if sline.startswith(section.headerTag):
              section.header = sline[1:].rstrip()
              continue

            if sline.startswith(">>END_MODULE"):
              section.inSection = False
              continue

            section.data[sampleName][datafilePrefix].append(sline.rstrip())
            continue

output_prefix = "error." + args.output if has_error else args.output

with open(output_prefix + ".summary.tsv", "w") as fout:
  fout.write("Sample\tFile\tCategory\tQCResult\n")
  for skey, svalue in sorted(summary.items()):
    slen = len(svalue)
    for vkey, vvalue in sorted(svalue.items()):
      for ckey, cvalue in sorted(vvalue.items()):
        fout.write("%s\t%s\t%s\t%s\n" % (skey, skey if slen==1 else vkey, ckey, cvalue))
    
with open(output_prefix + ".reads.tsv", "w") as fout:
  fout.write("Sample\tFile\tReads\n")
  for skey, svalue in sorted(reads.items()):
    slen = len(svalue)
    for vkey, vvalue in sorted(svalue.items()):
      fout.write("%s\t%s\t%s\n" % (skey, skey if slen==1 else vkey, vvalue))

with open(output_prefix + ".overrepresented.tsv", "w") as fout:
  fout.write("Sample\tiFile\t%s\n" % overrepresentedHeader)
  for skey, svalue in sorted(overrepresented.items()):
    slen = len(svalue)
    for vkey, vvalue in sorted(svalue.items()):
      fout.write("%s\t%s\t%s\n" % (skey, skey if slen==1 else vkey, vvalue))

for section in sections:
  with open("%s.%s.tsv" % (output_prefix, section.fileName), "w") as fout:
    fout.write("Sample\tFile\t%s\n" % section.header)
    for skey, svalue in sorted(section.data.items()):
      slen = len(svalue)
      for vkey, vvalue in sorted(svalue.items()):
        for avalue in vvalue:
          fout.write("%s\t%s\t%s\n" % (skey, skey if slen==1 else vkey, avalue))

