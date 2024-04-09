import pysam
import re
import logging
import argparse

from CountXmlUtils import readCountXmlQueriesInFeatures

DEBUG = False
NOT_DEBUG= not DEBUG

if DEBUG:
  input_xml="/nobackup/vickers_lab/projects/20230131_9074_ES_ARMseq_human_byMars_hg38/intermediate_data/bowtie1_genome_1mm_NTA_smallRNA_count/result/CAC_10_007DE_AlkB/CAC_10_007DE_AlkB.count.mapped.xml"
  input_bam = "/nobackup/vickers_lab/projects/20230608_9074_ES_ARMseq_human_byMars_HERV/intermediate_data/bowtie1_genome_1mm_NTA/result/CAC_10_007DE_AlkB.bam"
  name_pattern = "Arg-ACG"
  min_count = 2
  output_bam="/nobackup/vickers_lab/projects/20230131_9074_ES_ARMseq_human_byMars_hg38/bamsnap_locus/result/CAC_10_007DE_AlkB.filtered.bam"
else:
  parser = argparse.ArgumentParser(description="Generate smallRNA BAM from mapped xml.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input_xml', action='store', nargs='?', help='Input count xml file', required=NOT_DEBUG)
  parser.add_argument('-b', '--input_bam', action='store', nargs='?', help="Original bam file", required=NOT_DEBUG)
  parser.add_argument('-n', '--name_pattern', action='store', nargs='?', help="Name pattern", required=NOT_DEBUG)
  parser.add_argument('-c', '--min_count', action='store', nargs='?', type=int, default=2, help="Minimum read count of the query", required=NOT_DEBUG)
  parser.add_argument('-o', '--output_bam', action='store', nargs='?', help="Output bam file", required=NOT_DEBUG)

  args = parser.parse_args()
  
  print(args)
  
  input_xml = args.input_xml
  input_bam = args.input_bam
  name_pattern = args.name_pattern
  min_count = args.min_count
  output_bam = args.output_bam
  
logger = logging.getLogger('xmlToBamFilterName')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def accept_namepattern(feature_name):
  return re.search(name_pattern, feature_name)

logger.info("Reading smRNA from " + input_xml + " ...")
tDR_set = readCountXmlQueriesInFeatures(input_xml, accept_namepattern)
logger.info(f"There are {len(tDR_set)} queries matched with feature name pattern {name_pattern}")

qnames = [t.Name for t in tDR_set if t.Count >= min_count]
logger.info(f"There are {len(qnames)} queries with read count >= {min_count}")

logger.info(f"Saving reads to {output_bam} ...")
with pysam.Samfile(input_bam) as sam:
  header = sam.header

  with pysam.AlignmentFile(output_bam, "wb", header=header) as outf:
    for read in sam.fetch():
      if read.query_name in qnames:
        outf.write(read)

logger.info(f"Indexing {output_bam} ...")
pysam.index(output_bam)

logger.info(f"Done.")

