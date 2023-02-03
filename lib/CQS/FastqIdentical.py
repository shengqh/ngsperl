import logging
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

DEBUG = True

if DEBUG:
  input_file="/scratch/cqs/shengq2/ravi_shah_projects/20230201_apd_smallrna_hg38/preprocessing/identical_summary/result/apd_smallrna_hg38__fileList1.list"
  output_file="/scratch/cqs/shengq2/ravi_shah_projects/20230201_apd_smallrna_hg38/preprocessing/identical_summary/result/apd_smallrna_hg38.txt"
else:
  parser = argparse.ArgumentParser(description="Identical read summary.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input identical list file', required=True)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output summary file", required=True)

  args = parser.parse_args()
  
  print(args)
  
  input_file = args.input
  output_file = args.output

logger = logging.getLogger('identical')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("Reading identical table from " + input_file)

input_map = {}
with open(input_file, "rt") as fin:
  for line in fin:
    parts = line.rstrip().split('\t')
    input_map[parts[1]] = parts[0]

#print(input_map)

sample_names = sorted(input_map.keys())
n_reads = []
all_df = None
for sample_name in sample_names[0:3]:
  dup_file = input_map[sample_name]
  logger.info(f"reading {dup_file} ..." )

  #dup_file = "/scratch/cqs/ravi_shah_projects/shengq2/20230201_apd_smallrna_hg38/preprocessing/identical/result/CT_01_clipped_identical.fastq.dupcount"
  dup = pd.read_csv(dup_file, sep="\t")

  logger.info(f"  has {dup.shape[0]} unique reads" )
  n_reads.append(dup.shape[0])

  dup=dup.head(10000)
  dup['SeqLength'] = dup.Sequence.str.len()
  dup['logCount'] = np.log2(dup.Count)

  dup=dup[["SeqLength", "logCount"]]
  dup['Sample'] = sample_name

  if all_df is None:
    all_df = dup
  else:
    all_df = pd.concat([all_df, dup])

  sns.kdeplot(x = dup.SeqLength, y= dup.logCount, color='r', fill=True, cmap="Reds", thresh=0.05)  
  plt.savefig(sample_name + ".png", dpi=300)

n_read_df = pd.DataFrame({"Sample":sample_names, "UniqueRead":n_reads})
n_read_df.to_csv(output_file, sep="\t")

logger.info("Done")
