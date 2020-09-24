import argparse
import sys
import logging
import os
import subprocess
import gzip
import shutil

def gzip_file(logger, input_file, output_file):
  logger.info(f"compressing {input_file}")
  with open(input_file, 'rb') as f_in:
    with gzip.open(output_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
  logger.info(f"{output_file} saved.")
  return(output_file)

def run_cmd(logger, cmd):
  logger.info(cmd)
  subprocess.call(cmd, shell=True)

def convert(logger, i1, read1, read2, whitelist, output_prefix):
  res_r1_file = f"{output_prefix}_S1_L001_R1_001.fastq.gz"
  res_r2_file = f"{output_prefix}_S1_L001_R2_001.fastq.gz"

  cvt_r1_file_unzipped = "R1_sl.fastq.gz.4tenx.fastq"
  cvt_r2_file_unzipped = "R2_sl.fastq.gz.4tenx.fastq"

  cvt_r1_file = "R1_sl.fastq.gz.4tenx.fastq.gz"
  cvt_r2_file = "R2_sl.fastq.gz.4tenx.fastq.gz"

  if not os.path.isfile(cvt_r1_file):
    if not os.path.isfile(cvt_r1_file_unzipped):
      cmd = f"ust10x -wl {whitelist} -i1 {i1} -r1 {read1} -r2 {read2}"
      run_cmd(logger, cmd)
    gzip_file(logger, cvt_r1_file_unzipped, cvt_r1_file)
    gzip_file(logger, cvt_r2_file_unzipped, cvt_r2_file)
    os.delete(cvt_r1_file_unzipped)
    os.delete(cvt_r2_file_unzipped)
  
  os.rename(cvt_r1_file_unzipped, res_r1_file)
  os.rename(cvt_r2_file_unzipped, res_r2_file)

  logger.info("done")
    
def main():
  DEBUG=False
  NotDEBUG=not DEBUG

  parser = argparse.ArgumentParser(description="Perform ust10x and rename rsult",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i1', '--i1', action='store', nargs='?', help='Input I1 file', required=NotDEBUG)
  parser.add_argument('-r1', '--read1', action='store', nargs='?', help='Input R1 file', required=NotDEBUG)
  parser.add_argument('-r2', '--read2', action='store', nargs='?', help='Input R2 file', required=NotDEBUG)
  parser.add_argument('-wl', '--whitelist', action='store', nargs='?', help='Input R2 file', required=NotDEBUG)
  parser.add_argument('-o', '--output_prefix', action='store', nargs='?', help="Output file prefix", required=NotDEBUG)

  args = parser.parse_args()

  if DEBUG:
    args.i1 = "/gpfs52/data/vickers_lab/20200828_5059_ES_tellseq/tellseq_read/Nova239TellSeq_I1_T501.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz"
    args.read1 = "/gpfs52/data/vickers_lab/20200828_5059_ES_tellseq/tellseq_read/Nova239TellSeq_R1_T501.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz"
    args.read2 = "/gpfs52/data/vickers_lab/20200828_5059_ES_tellseq/tellseq_read/Nova239TellSeq_R2_T501.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz"
    args.whitelist = "/scratch/cqs_share/softwares/tellseq/conversion_tool/4M-with-alts-february-2016.txt"
    args.output_prefix = "T501"

  logger = logging.getLogger('ust10x')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  convert(logger, args.i1, args.read1, args.read2, args.whitelist, args.output_prefix)

if __name__ == "__main__":
    main()
