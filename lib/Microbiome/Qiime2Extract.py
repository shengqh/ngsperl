import os
import argparse
import logging
from asyncore import read
from zipfile import ZipFile
from Qiime2Utils import extract_file 

def main():
  parser = argparse.ArgumentParser(description="Extract file from qiime2 qzv file.",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Qiime2 qzv file', required=True)
  parser.add_argument('-p', '--path', action='store', nargs='?', help="Input path in zipped file", required=True)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file name")

  args = parser.parse_args()

  logger = logging.getLogger('extract_file')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

  extract_file(args.input, args.path, args.output, logger)

if __name__ == "__main__":
  main()

