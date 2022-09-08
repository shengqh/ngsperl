import argparse
import logging
import os

#https://software.broadinstitute.org/gatk/documentation/article.php?id=6926

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Prepare multi fastq file list for GATK4",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', required=NotDEBUG, help='Input FASTQ file')
parser.add_argument('-n', '--name', action='store', nargs='?', required=NotDEBUG, help='Input sample name')
parser.add_argument('-j', '--json', action='store', nargs='?', required=NotDEBUG, help='Input json template file')
parser.add_argument('-o', '--output', action='store', nargs='?', required=NotDEBUG, help='Output list file')

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-581_1.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-581_2.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-582_1.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-582_2.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-583_1.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-583_2.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-584_1.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-584_2.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-585_1.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-585_2.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-586_1.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-586_2.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-587_1.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-587_2.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-588_1.fq.gz,/scratch/cqs/breast_cancer_spore/clean/01-127/V350089367_L04_B5GHUMtkywRAADA-588_2.fq.gz"
  args.name = "P01_1271"
  args.json = "/data/cqs/softwares/cqsperl/config/wdl/processing-for-variant-discovery-gatk4.hg38.wgs.inputs.json"
  args.output = "/scratch/cqs/breast_cancer_spore/analysis/prepare_config/result/P01_127/P01_127"

print(args)

logger = logging.getLogger('prepare')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

files = args.input.split(',')

flow_list_file = args.output + ".txt"
with open(flow_list_file, "wt") as flist:
  for idx in range(0, len(files), 2):
    fastq_list_file = f"{args.output}_{str(idx).zfill(3)}.list"
    flist.write(os.path.abspath(fastq_list_file) + "\n")
    with open(fastq_list_file, "wt") as fout:
      fout.write(files[idx] + "\n")
      fout.write(files[idx+1] + "\n")

def replace_value(line, newvalue):
  parts = line.split(":", 1)
  parts[1] = f" \"{newvalue}\",\n"
  return(":".join(parts))

json_file = args.output + ".inputs.json"
with open(json_file, "wt") as fjson:
  with open(args.json, "rt") as fin:
    for line in fin:
      if "PreProcessingForVariantDiscovery_GATK4.sample_name" in line:
        line = replace_value(line, args.name)
      elif "PreProcessingForVariantDiscovery_GATK4.flowcell_unmapped_bams_list" in line:
        line = replace_value(line, flow_list_file)
      elif "PreProcessingForVariantDiscovery_GATK4.unmapped_bam_suffix" in line:
        line = replace_value(line, ".list")
      fjson.write(line)

logger.info("done.")
