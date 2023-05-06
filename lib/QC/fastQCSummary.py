import argparse
import sys
import logging
import os

from collections import defaultdict

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

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Summarize the fastqc result.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input fastqc data file list', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output summary prefix", required=NotDEBUG)

args = parser.parse_args()

if DEBUG:
  args.input = "/nobackup/vickers_lab/projects/20230502_9880_smallRNA_rice_hg38_byTiger/preprocessing/fastqc_raw_summary/result/fileList1.txt"
  args.output = "/nobackup/vickers_lab/projects/20230502_9880_smallRNA_rice_hg38_byTiger/preprocessing/fastqc_raw_summary/result/P9880.FastQC"

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

sections = [Section("baseQuality", ">>Per base sequence quality", "#Base"),
            Section("baseContent", ">>Per base sequence content", "#Base"),
            Section("sequenceGC", ">>Per sequence GC content", "#GC Content")]
summary={}
reads={}
overrepresented={}
overrepresentedHeader = ""
adapter=defaultdict(dict)
all_adapters = {}

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
      logger.info(f"reading {datafile} ...")
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


        if sline.startswith(">>Adapter Content"):
          bInAdapter = True
          continue

        if bInAdapter:
          if sline.startswith("#Position"):
            adapter_headers = sline[1:].rstrip().split('\t')
            adapter[sampleName][datafilePrefix] = []
            for ah in adapter_headers[1:]:
              all_adapters[ah] = 1
            continue
          if not sline.startswith(">>END_MODULE"):
            parts = sline.rstrip().split('\t')
            ad_dic = {}
            for idx in range(1, len(adapter_headers)):
              ad_dic[adapter_headers[idx]] = parts[idx]
            adapter[sampleName][datafilePrefix].append([parts[0], ad_dic])
          else:
            bInOver = False
          continue

        for section in sections:
          if sline.startswith(section.sectionTag):
            section.inSection = True
            section.data[sampleName][datafilePrefix] = []
            continue
          
          if section.inSection:
            if sline.startswith(section.headerTag):
              cur_header = sline[1:].rstrip()
              if section.header == "":
                section.header = cur_header
              elif section.header != cur_header:
                raise Exception(f"Header not match:\n{section.header}\n{cur_header}\n")
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
  fout.write("Sample\tFile\t%s\n" % overrepresentedHeader)
  for skey, svalue in sorted(overrepresented.items()):
    slen = len(svalue)
    for vkey, vvalue in sorted(svalue.items()):
      fout.write("%s\t%s\t%s\n" % (skey, skey if slen==1 else vkey, vvalue))

with open(output_prefix + ".adapter.tsv", "w") as fout:
  all_adapters_names = sorted(all_adapters.keys())
  all_adapters_names = [an for an in all_adapters_names if not "Poly" in an]
  fout.write("Sample\tFile\tPosition\t%s\n" % "\t".join(all_adapters_names))
  for skey, svalue in sorted(adapter.items()):
    slen = len(svalue)
    for vkey, vvalue in sorted(svalue.items()):
      sample = skey if slen==1 else vkey
      for pvalue in vvalue:
        pos = pvalue[0]
        pos_dic = pvalue[1]
        fout.write(f"{skey}\t{sample}\t{pos}")
        for an in all_adapters_names:
          if an in pos_dic:
            fout.write(f"\t{pos_dic[an]}")
          else:
            fout.write(f"\t0.0")
        fout.write("\n")

for section in sections:
  with open("%s.%s.tsv" % (output_prefix, section.fileName), "w") as fout:
    fout.write("Sample\tFile\t%s\n" % section.header)
    for skey, svalue in sorted(section.data.items()):
      slen = len(svalue)
      for vkey, vvalue in sorted(svalue.items()):
        for avalue in vvalue:
          fout.write("%s\t%s\t%s\n" % (skey, skey if slen==1 else vkey, avalue))

