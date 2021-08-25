import logging
import argparse
import subprocess
import os
from FileListUtils import readFileMap

DEBUG = 0

if DEBUG:
  inputFile="/scratch/vickers_lab/projects/20210504_6282_RA_smRNAseq_mouse_byTiger/covid19/bowtie1_covid19_04_table/result/RA_6282_mouse__fileList1.list"
  outputFile="/scratch/vickers_lab/projects/20210504_6282_RA_smRNAseq_mouse_byTiger/covid19/bowtie1_covid19_04_table/result/RA_6282_mouse.count"
  name_map_file=None
else:
  parser = argparse.ArgumentParser(description="Build count table.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count list file')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output result file")

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  outputFile = args.output

logger = logging.getLogger('count')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("Reading count table from " + inputFile)
countList = readFileMap(inputFile)

sample_names = sorted(countList.keys())

count_map = {}
chroms = set()
for sample_name in sample_names:
  sample_file = countList[sample_name]
  file_map = {}
  count_map[sample_name] = file_map
  with open(sample_file, 'rt') as fin:
    fin.readline()
    for line in fin:
      parts = line.rstrip().split('\t')
      file_map[parts[0]] = [parts[1], parts[2]]
      chroms.add(parts[0])

sorted_chroms = sorted(chroms)
count_names=["unique", "read"]

with open(outputFile, "wt") as fout:
  fout.write("Feature\t%s\n" % "\t".join(sample_names))
  for chrom in sorted_chroms:
    for idx in [0,1]:
      fout.write(chrom + "_" + count_names[idx])
      for sample_name in sample_names:
        file_map = count_map[sample_name]
        if chrom in file_map:
          fout.write("\t%s" % file_map[chrom][idx])
        else:
          fout.write("\t0")
      fout.write("\n")

rscript = os.path.realpath(__file__) + ".R"
subprocess.call ("R --vanilla -f " + rscript + " --args \"" + outputFile + "\" \"" + outputFile + ".png\"", shell=True)

logger.info("Done")
