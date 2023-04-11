import argparse
import logging
import gzip
import os

def validate2(logger, input1str, input2str, output):
  input1 = input1str.split(',')
  input2 = input2str.split(',')

  error_files = []
  error_msgs = []
  total_read_count = 0
  idx = 0
  while idx < len(input1):
    read1 = input1[idx]
    read2 = input2[idx]
    idx += 1

    logger.info("Validating %s ..." % read1)
    logger.info("Validating %s ..." % read2)
    read_count = 0
    error_msg = None
    try:
      if read1.endswith(".gz"):
        fin1 = gzip.open(read1, "rt")
        fin2 = gzip.open(read2, "rt")
      else:
        fin1 = open(read1, "rt")
        fin2 = open(read2, "rt")

      with fin1, fin2:
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
            #break

          line1 = line1.rstrip()
          line2 = line2.rstrip()

          qname1 = line1.split(' ', 1)[0]
          qname2 = line2.split(' ', 1)[0]
          if qname1 != qname2:
            if qname1.endswith("/1"):
              qname1 = qname1[0:len(qname1) - 2]
              qname2 = qname2[0:len(qname2) - 2]
              if qname1 != qname2:
                error_msg = "query name not equals: %s , %s in %s and %s" % (qname1, qname2, read1, read2)
                break
            else:
              error_msg = "query name not equals: %s , %s in %s and %s" % (qname1, qname2, read1, read2)
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
    except OSError as e:
      error_msg = "I/O error({0}): {1}".format(e.errno, e.strerror)
    except Exception as e:
      error_msg = "Unexpected error: {0}".format(e)

    total_read_count += read_count

    if error_msg != None:
      logger.error(error_msg)
      error_msgs.append(error_msg)
      error_files.append(read1)
      error_files.append(read2)

  error_file = output + ".error"
  if os.path.exists(error_file):
    os.remove(error_file)

  if len(error_msgs) > 0:
    logger.error("failed : all error msgs:")
    if os.path.exists(output):
      os.remove(output)
    with open(error_file, "wt") as fout:
      for error_msg in error_msgs:
        logger.error(error_msg)
        fout.write("ERROR: %s\n" % error_msg)

      fout.write("\n")
      for error_file in error_files:
        fout.write(error_file + "\n")
      
    return(1)
  else:
    logger.info("succeed")
    with open(output, "wt") as fout:
      fout.write("READ\t%d\n" % total_read_count)
    return(0)

def validate1(logger, input, output): 
  parts = input.split(',')
  input1 = []
  input2 = []
  idx = 0
  while idx < len(parts):
    input1.append(parts[idx])
    input2.append(parts[idx+1])
    idx += 2
  
  return(validate2(logger, ",".join(input1), ",".join(input2), output))

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Validate pairend FASTQ.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input fastq files, joined by ","', required=NOT_DEBUG)
  parser.add_argument('--input2', action='store', nargs='?', help='Input fastq read2 files if input contains read1 only, joined by ","')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NOT_DEBUG)
  
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/cqs/shengq2/vivian_weiss/20200831_human_wntpathway_rnaseq/cutadapt/result/X5909CP3_clipped.1.fastq.gz,/scratch/cqs/shengq2/vivian_weiss/20200831_human_wntpathway_rnaseq/cutadapt/result/X5909CP3_clipped.2.fastq.gz"
    args.output = "/scratch/cqs/shengq2/vivian_weiss/20200831_human_wntpathway_rnaseq/cutadapt/X5909CP3.txt"
  
  logger = logging.getLogger('fastq_validator')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  if args.input2 is None:
    return(validate1(logger, args.input, args.output))
  else:
    return(validate2(logger, args.input, args.input2, args.output))
  
if __name__ == "__main__":
  main()
