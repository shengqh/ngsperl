import argparse
import sys
import logging
import os
import csv
import warnings

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Post annovar to maf processor.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input MAF file for single sample', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output MAF file for multiple sample", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/cqs/shengq2/macrae_linton/20180913_linton_exomeseq_2118_human_cutadapt/bwa_refine_gatk4_hc_gvcf_vqsr_filterMAF_annovar_filter_toMAF/result/linton_exomeseq_2118.freq0.001.filtered.tsv.tmp"
  args.output = "/scratch/cqs/shengq2/macrae_linton/20180913_linton_exomeseq_2118_human_cutadapt/bwa_refine_gatk4_hc_gvcf_vqsr_filterMAF_annovar_filter_toMAF/result/linton_exomeseq_2118.freq0.001.filtered.tsv.maf"

logger = logging.getLogger('annovar2maf')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

with open(args.input, "r") as fin:
  mafheaders = fin.readline().rstrip().split("\t")
  sample_index = mafheaders.index("FORMAT") + 1
  Tumor_Sample_Barcode_index = mafheaders.index("Tumor_Sample_Barcode")
  Reference_Allele_index = mafheaders.index("Reference_Allele")
  Tumor_Seq_Allele2_index = mafheaders.index("Tumor_Seq_Allele2")
  
  tmpFile = args.output + ".tmp"
  newheaders = mafheaders[0:(sample_index-1)]
  with open(tmpFile, "w") as fout:
    fout.write('\t'.join(newheaders) + "\tsample_id\tt_depth\tt_ref_count\tt_alt_count\tFORMAT\tGENOTYPE\n")
    for line in fin:
      maf_items = line.rstrip().split("\t")
      new_items = maf_items[0: (sample_index-1)]

      #get AD and DP column
      sampleFormats=maf_items[(sample_index-1)]
      sampleFormatsSplit = sampleFormats.split(":")
      if "DP" in sampleFormatsSplit:
        format_DP_index=sampleFormatsSplit.index("DP")
      else:
        format_DP_index=2
      if "AD" in sampleFormatsSplit:
        format_AD_index=sampleFormatsSplit.index("AD")
      else:
        format_AD_index="NA"
        #warnings.warn("Can't find AD in " + line)
        #continue
      
      for sample_idx in range(sample_index, len(mafheaders)):
        sample_name = mafheaders[sample_idx]
        new_items[Tumor_Sample_Barcode_index] = sample_name

        sample_data = maf_items[sample_idx]

        if sample_data.startswith("0/0") or sample_data.startswith("0|0") or sample_data.startswith("./.")  or sample_data.startswith(".|.") :
          continue

        if not (sample_data.startswith("0/1:") or sample_data.startswith("0|1:") or sample_data.startswith("1/0:") or sample_data.startswith("1|0:") or sample_data.startswith("1/1:") or sample_data.startswith("1|1:")):
          raise Exception('I don\'t know genotype: ' + sample_data)
        
        parts = sample_data.split(":")
        t_depth = parts[format_DP_index]
        if (format_AD_index=="NA"):
          t_af = float(parts[format_AF_index])
          t_alt_count=str(int(int(t_depth)*float(t_af)))
          t_ref_count=str(int(t_depth)-int(t_alt_count))
        else:
          alleles = parts[format_AD_index].split(",")
          t_ref_count = alleles[0]
          t_alt_count = alleles[1]

        fout.write("\t".join(new_items) + "\t" + sample_name + "\t" + t_depth + "\t" + t_ref_count + "\t" + t_alt_count + "\t" + maf_items[sample_index-1] + "\t" + sample_data + "\n")
    
  if os.path.isfile(args.output):
    os.remove(args.output)
  os.rename(tmpFile, args.output)
