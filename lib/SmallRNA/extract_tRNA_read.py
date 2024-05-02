import logging
import argparse
import os
import gzip
import xml.etree.ElementTree as ET
from FileListUtils import readFileMap

def accept_func(feature_name:str):
  return(feature_name.startswith("tRNA:"))

# assert(not accept_func("miRNA:hsa-let-7a-5p"))
# assert(accept_func("tRNA:hsa-let-7a-5p"))

def extract(logger, input_list_file, output_file, add_cca):
  xml_files = readFileMap(input_list_file)
  sample_names = sorted(xml_files.keys())

  with gzip.open(output_file, "wt") as fout:
    fout.write("sample\tcategory\tfeature_name\tquery_name\tmapped_sequence\tmapped_length\toriginal_sequence\toriginal_length\tmap_offset\tquery_count\tnum_mismatch\n")

    icount = 0
    for sample_name in sample_names:
      icount += 1
      # if icount % 5 == 0:
      #   break

      xml_file = xml_files[sample_name]

      logger.info(f"Processing {xml_file} ...")

      tree = ET.parse(xml_file)
      root = tree.getroot()
      query_seq_map = {}

      queries_node = root.find('queries')
      for query_node in queries_node.findall('query'):
        query_name = query_node.get("name").rstrip()
        if "sequence" in query_node.attrib:
          query_sequence=query_node.get("sequence")
        else:
          query_sequence=query_node.get("seq")
        query_seq_map[query_name] = query_sequence

      #print(query_seq_map)
      result_node = root.find('subjectResult')
      feature_map = {}

      for feature_group_node in result_node.findall('subjectGroup'):
        feature_nodes = feature_group_node.findall('subject')
        for feature_node in feature_nodes:
          feature_name = feature_node.get("name")
          if not accept_func(feature_name):
            continue
          
          region_nodes = feature_node.findall('region')
          for region_node in region_nodes:
            seqname = region_node.get("seqname")
            if not seqname.startswith('Homo_sapiens_'):
              continue

            sequence = region_node.get("sequence")
            cur_queries = {}

            query_nodes = region_node.findall('query')
            for query_node in query_nodes:
              qname = query_node.get('qname')
              qseq = query_seq_map[qname]
              cur_queries[qname] = {
                "sample": sample_name,
                "offset": query_node.get('offset'),
                "query_count": query_node.get('query_count'),
                "num_mismatch": query_node.get('nmi'),
                "sequence": qseq
              }

            feature_map[feature_name] = {
              "sequence": sequence,
              "queries": cur_queries 
            }

      for feature_name in sorted(feature_map.keys()):
        fmap = feature_map[feature_name]

        short_feature_name = feature_name.replace("tRNA:","")
        fseq_cca = fmap['sequence']
        if add_cca:
          fseq_cca = fseq_cca + "CCA"
        fout.write(f"{sample_name}\tparent\t{short_feature_name}\t\t{fmap['sequence']}\t{len(fmap['sequence'])}\t{fseq_cca}\t{len(fseq_cca)}\t0\t0\t0\n")
        for query_name in fmap['queries'].keys():
          qmap = fmap['queries'][query_name]
          trimmed_sequence = query_name.split("CLIP_")[1]
          original_sequence = qmap['sequence'] + trimmed_sequence
          fout.write(f"{sample_name}\tread\t{short_feature_name}\t{query_name}\t{qmap['sequence']}\t{len(qmap['sequence'])}\t{original_sequence}\t{len(original_sequence)}\t{qmap['offset']}\t{qmap['query_count']}\t{qmap['num_mismatch']}\n")

DEBUG = False
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="Extract tRNA read from Xml",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count xml list file', required=NOT_DEBUG)
parser.add_argument('--add_cca', action='store_true', help='Add CCA to parent sequence')
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output tRNA file", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="/nobackup/vickers_lab/projects/20221122_9074_ES_ARMseq_human_byMars/intermediate_data/host_tRNA_mismatch_vis/result/9074_ES__fileList1.list"
  args.add_cca=True
  args.output="/nobackup/vickers_lab/projects/20221122_9074_ES_ARMseq_human_byMars/intermediate_data/host_tRNA_mismatch_vis/result/9074_ES.tRNA.txt"
print(args)
  
logger = logging.getLogger('extract_tRNA')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

if not os.path.exists(args.output):
  extract(logger, args.input, args.output, args.add_cca)

logger.info("done.")
