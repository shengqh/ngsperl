import argparse
import logging
import gzip

def validate(logger, input, output): 
  parts = input.split(',')
  read1 = parts[0]
  read2 = parts[1]
  logger.info("Validating %s ..." % read1)
  logger.info("Validating %s ..." % read2)
  read_count = 0
  error_msg = ""
  with gzip.open(read1, "rt") as fin1:
    with gzip.open(read2, "rt") as fin2:
      while(True):
        line1 = fin1.readline()
        line2 = fin2.readline()
        if not line1:
          if not line2:
            break
          else:
            error_msg = "%s end but %s not" % (read1, read2)
            break
        if not line2:
          error_msg = "%s end but %s not" % (read2, read1)
          break

        if not line1.startswith('@'):
          error_msg = "query %s not starts with @ in %s" % (line1, read1)
          break

        if not line2.startswith('@'):
          error_msg = "query %s not starts with @ in %s" % (line2, read2)
          break

        read_count += 1
        if read_count % 1000000 == 0:
          logger.info("%s" % read_count)

        qname1 = line1.split(' ', 1)[0]
        qname2 = line2.split(' ', 1)[0]
        if qname1 != qname2:
          error_msg = "query name not equals: %s , %s" % (qname1, qname2)
          break

        seq1 = fin1.readline()
        mid1 = fin1.readline()
        score1 = fin1.readline()

        if not seq1 or not mid1 or not score1:
          error_msg = "unexpected end for query %s in %s" % (qname1, read1)
          break

        if len(seq1) != len(score1):
          error_msg = "sequence length not equals to score length for query %s in %s\n  seq  :%s\n  score:%s\n" % (qname1, read1, seq1, score1)
          break

        seq2 = fin2.readline()
        mid2 = fin2.readline()
        score2 = fin2.readline()

        if not seq2 or not mid2 or not score2:
          error_msg = "unexpected end for query %s in %s" % (qname2, read2)
          break

        if len(seq2) != len(score2):
          error_msg = "sequence length not equals to score length for query %s in %s\n  seq  :%s\n  score:%s\n" % (qname2, read2, seq2, score2)
          break
    
    if error_msg != "":
      with open(output + ".error", "wt") as fout:
        fout.write("ERROR: %s" % error_msg)
      logger.error(error_msg)
      return(1)
    else:
      with open(output, "wt") as fout:
        fout.write("READ\t%d" % read_count)
      logger.info("done")
      return(0)

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Validate pairend FASTQ.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input fastq files, joined by ","', required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NOT_DEBUG)
  
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/cqs/shengq2/vivian_weiss/20200831_human_wntpathway_rnaseq/cutadapt/result/X5909CP3_clipped.1.fastq.gz,/scratch/cqs/shengq2/vivian_weiss/20200831_human_wntpathway_rnaseq/cutadapt/result/X5909CP3_clipped.2.fastq.gz"
    args.output = "/scratch/cqs/shengq2/vivian_weiss/20200831_human_wntpathway_rnaseq/cutadapt/X5909CP3.txt"
  
  logger = logging.getLogger('fastq_validator')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  return(validate(logger, args.input, args.output))
  
if __name__ == "__main__":
  main()
