import logging
import argparse
import subprocess
import os
import pickle
import statistics
from FileListUtils import readFileMap

DEBUG = 0

if DEBUG:
  inputFile="/scratch/vickers_lab/projects/20210817_smallRNA_contamination/3626_RA/bowtie1_contamination_all_03_table/result/3626_RA__fileList1.list"
  objectFile='/scratch/vickers_lab/projects/20210817_smallRNA_contamination/3626_RA/bowtie1_contamination_all_03_table/result/3626_RA__fileList2.list'
  outputFile="/scratch/vickers_lab/projects/20210817_smallRNA_contamination/3626_RA/bowtie1_contamination_all_03_table/result/3626_RA.count.txt"
else:
  parser = argparse.ArgumentParser(description="Build count table.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count list file', required=True)
  parser.add_argument('-s', '--object', action='store', nargs='?', help='Input count object list file', required=True)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output result file", required=True)

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  objectFile = args.object
  outputFile = args.output

logger = logging.getLogger('count')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')


logger.info("Reading count map table from " + objectFile)
objectList = readFileMap(objectFile)

sample_names = sorted(objectList.keys())
sample_objs = {}
for sample in sample_names:
  object_file = objectList[sample]
  with (open(object_file, "rb")) as fin:
    sample_objs[sample] = pickle.load(fin)

chroms = set()
seqcount={}
for sample in sample_names:
  count_map = sample_objs[sample]
  seq_map = {}
  for chrom in count_map.keys():
    chroms.add(chrom)
    chrom_map = count_map[chrom]
    for seq in chrom_map.keys():
      seqcount.setdefault(seq, [{}, {sample:0 for sample in sample_names}])
      seqcount[seq][0][chrom] = f"{chrom_map[seq][1]}-{chrom_map[seq][2]}"
      seqcount[seq][1][sample] = chrom_map[seq][0]

#calculate seq median value cross samples
for seq in seqcount.keys():
  median_count = statistics.median(seqcount[seq][1].values())
  seqcount[seq].append(median_count)

sorted_seqcount = sorted(seqcount.items(), key=lambda item: item[1][2],  reverse=True)
sorted_chroms = sorted(chroms)

count = 0
with open(outputFile + ".top50.txt", "wt") as fout:
  fout.write("Sequence\t%s\t%s\n" % ("\t".join(sorted_chroms), "\t".join(sample_names)))
  for item in sorted_seqcount:
    seq = item[0]
    chroms = item[1][0]
    count_map = item[1][1]
    fout.write(seq)

    for chrom in sorted_chroms:
      if chrom in chroms:
        fout.write(f"\t{chrom}:{chroms[chrom]}")
      else:
        fout.write("\t")

    for sample in sample_names:
      fout.write("\t%d" % count_map[sample])
    
    fout.write("\n")

    count += 1
    if count == 50:
      break

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
