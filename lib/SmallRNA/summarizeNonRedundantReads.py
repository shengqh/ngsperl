import sys
import os
import logging
import argparse
import xml.etree.ElementTree as ET
from CountXmlUtils import readCountXmlQueriesInFeatures
from DupCountUtils import readDupCountQueries
from os.path import basename

DEBUG = True
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="Summarize Non-Redundant Reads",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count xml file', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output bam file", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="W:/SequenceProjects/3018-KCV-77_78_79/20191014_smallRNA_3018-KCV-77_78_79_mouse_v4_liverWT_byTiger/host_genome/bowtie1_genome_1mm_NTA_smallRNA_count/result/Liver_WT_12/Liver_WT_12.count.mapped.xml"
  args.nodup="W:/SequenceProjects/3018-KCV-77_78_79/20191014_smallRNA_3018-KCV-77_78_79_mouse_v4_liverWT_byTiger/preprocessing/identical/result/Liver_WT_12_clipped_identical.fastq.dupcount"
  args.output="W:/SequenceProjects/3018-KCV-77_78_79/20191014_smallRNA_3018-KCV-77_78_79_mouse_v4_liverWT_byTiger/data_visualization/nonredundant_summary/nonredundant_summary.txt"
  
print(args)
  
logger = logging.getLogger('summarizeNonRedundantReads')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

sampleNames = ["Liver_WT_01", "Liver_WT_03", "Liver_WT_05", "Liver_WT_08", "Liver_WT_10", "Liver_WT_12", "Liver_WT_14"]
  
def acceptMirna(featureName):
  return featureName.startswith("miRNA")

with open(args.output, "wt") as fout:
  fout.write("Sample\tTotalReads\tTotalUniqueReads\tTotalMirnaReads\tTotalUniqueMirnaReads\n")
  for sampleName in sampleNames:
    nodupfile = args.nodup.replace('Liver_WT_12', sampleName)
    inputfile = args.input.replace('Liver_WT_12', sampleName)
    
    logger.info("Reading " + nodupfile + " ...")
    queries = readDupCountQueries(nodupfile, 0, sorted=False)
  
    logger.info("Reading miRNA from " + inputfile + " ...")
    mirnas = readCountXmlQueriesInFeatures(inputfile, acceptMirna)

    fout.write("%s\t%d\t%d\t%d\t%d\n"%(sampleName, sum([q.Count for q in queries]), len(queries), sum([q.Count for q in mirnas]), len(mirnas) ))

logger.info("done.")
    
  
