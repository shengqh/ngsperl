import pysam
import argparse
import logging
import matplotlib.pyplot as plt

DEBUG = False
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="Get bam template length distribution",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file, sorted by query name', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output template length distribution file", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="/nobackup/brown_lab/projects/20241022_12329_cutrun_mm10_histone/bowtie2/result/IgG/IgG.rmdup.clean.bam"
  args.output="IgG.rmdup.clean.tlen.txt"

logger = logging.getLogger('bamTLEN')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def is_primary_read(read):
    return not read.is_secondary and not read.is_supplementary

length_dic = {}

with pysam.Samfile(args.input, "rb", check_sq=False) as sam:
  processed = 0

  for read in sam.fetch(until_eof=True):
    processed += 1
    if processed % 1000000 == 0:
      logger.info(f"{processed} reads processed.")

    if not is_primary_read(read):
      continue

    if read.is_reverse:
      continue
    
    #only count the primary and forward reads
    template_length = abs(read.tlen)
    if template_length not in length_dic:
      length_dic[template_length] = 1
    else:
      length_dic[template_length] += 1

sorted_lens = sorted(length_dic.keys())
sorted_values = [length_dic[l] for l in sorted_lens]

with open(args.output, "wt") as fout:
  fout.write("TLEN\tCount\n")
  for tlen in sorted_lens:
    fout.write(f"{tlen}\t{length_dic[tlen]}\n")

plt.rcParams["figure.figsize"] = (7, 4)
plt.bar(sorted_lens, sorted_values)
plt.xlabel('Template length of primary mapped reads')
plt.ylabel('Read count')
plt.savefig(args.output + '.png', dpi=300)

logger.info("done")