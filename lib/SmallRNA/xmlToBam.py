import pysam
import sys
import gzip
import os
import logging
import argparse
import xml.etree.ElementTree as ET

DEBUG = 0
NTA_TAG = ":CLIP_"

if DEBUG:
  inputFile="Z:/Shared/Labs/Vickers Lab/Tiger/projects/20161223_smallRNA_3018-KCV-77_78_79_86_v3.T/host_genome/bowtie1_genome_1mm_NTA_smallRNA_count/result/APOB_SRBIKO_86_01/APOB_SRBIKO_86_01.count.mapped.xml"
  outputFile="H:/temp/APOB_SRBIKO_86_01.bam"
else:
  parser = argparse.ArgumentParser(description="Generate smallRNA BAM from mapped xml.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Fastq file')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output NTA Fastq file")
  parser.set_defaults(ccaa=False)

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  outputFile = args.output

logger = logging.getLogger('fastqSmallRnaNTA')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

tree = ET.parse(inputFile)
root = tree.getroot()

queries = root.find('queries')
for query in queries.findall('query'):
  print(query.attrib)
# 
# with pysam.AlignmentFile(tmpfilename, "wb", header=header) as outf:
#     a = pysam.AlignedSegment()
#     a.query_name = "read_28833_29006_6945"
#     a.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
#     a.flag = 99
#     a.reference_id = 0
#     a.reference_start = 32
#     a.mapping_quality = 20
#     a.cigar = ((0,10), (2,1), (0,25))
#     a.next_reference_id = 0
#     a.next_reference_start=199
#     a.template_length=167
#     a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
#     a.tags = (("NM", 1),
#               ("RG", "L1"))
#     outf.write(a)
    