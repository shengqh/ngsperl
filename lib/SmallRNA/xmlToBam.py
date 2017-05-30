import pysam
import sys
import gzip
import os
import logging
import argparse
import xml.etree.ElementTree as ET
from Bio.Seq import Seq

DEBUG = 0

if DEBUG:
  inputFile="/scratch/cqs/shengq1/vickers/20161223_smallRNA_3018-KCV-77_78_79_86_v3.T/host_genome/bowtie1_genome_1mm_NTA_smallRNA_count/result/APOB_SRBIKO_86_01/APOB_SRBIKO_86_01.count.mapped.xml"
  oldbamFile = "/scratch/cqs/shengq1/vickers/20161223_smallRNA_3018-KCV-77_78_79_86_v3.T/host_genome/bowtie1_genome_1mm_NTA/result/APOB_SRBIKO_86_01/APOB_SRBIKO_86_01.bam"
  outputFilePrefix="/scratch/cqs/shengq1/temp/APOB_SRBIKO_86_01"
else:
  parser = argparse.ArgumentParser(description="Generate smallRNA BAM from mapped xml.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count xml file')
  parser.add_argument('-b', '--oldbamfile', action='store', nargs='?', help="Original bam file")
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output bam file")

  args = parser.parse_args()
  
  print(args)
  
  inputFile = args.input
  oldbamFile = args.oldbamfile
  outputFilePrefix = args.output

logger = logging.getLogger('xmlToBam')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

sam = pysam.Samfile(oldbamFile)
header = sam.header
chrMap = {}
for idx, sq in enumerate(header["SQ"]):
  chr = sq['SN']
  chrMap[chr] = idx
  if chr.startswith("chr"):
    chrMap[chr[3:]] = idx
  
tree = ET.parse(inputFile)
root = tree.getroot()

queries = root.find('queries')

unsorted = outputFilePrefix + ".unsorted.bam"
with pysam.AlignmentFile(unsorted, "wb", header=header) as outf:
  for query in queries.findall('query'):
    query_count = int(query.get("count")) + 1
    query_name = query.get("name")
    query_sequence=query.get("sequence")
    reverse_query_sequence=str(Seq(query_sequence).reverse_complement())
    #query_qualities=pysam.qualitystring_to_array("<" * len(query_sequence))
    for loc in query.findall("location"):
      strand = loc.get("strand")

      a = pysam.AlignedSegment()
      
      if(strand == '+'):
        a.flag = 0
        a.query_sequence=query_sequence
      else:
        a.flag = 16
        a.query_sequence=reverse_query_sequence
        
      a.reference_id = chrMap[ loc.get("seqname")]
      a.reference_start = int(loc.get("start")) - 1
      a.mapping_quality = 255
      a.cigarstring = loc.get("cigar")
      a.next_reference_id = -1
      a.next_reference_start=-1
      a.template_length=0
      #a.query_qualities = query_qualities
      a.tags = (("XA", int(loc.get("score"))),
                ("MD", loc.get("mdz")),
                ("NM", int(loc.get("nmi"))))

      for idx in range(1, query_count):
        a.query_name = query_name + ":" + str(idx)
        outf.write(a)
   
pysam.sort(unsorted, outputFilePrefix)
pysam.index(outputFilePrefix + ".bam")
os.remove(unsorted)