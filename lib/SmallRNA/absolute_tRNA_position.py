import sys
import os
import logging
import argparse
import time
import xml.etree.ElementTree as ET

# xml_file="/nobackup/vickers_lab/projects/20221122_9074_ES_ARMseq_human_byMars/intermediate_data/bowtie1_genome_1mm_NTA_smallRNA_count/result/CAC_10_007DE/CAC_10_007DE.count.mapped.xml"
# position_file="/nobackup/vickers_lab/projects/20221122_9074_ES_ARMseq_human_byMars_tRNA_pos/trna_absolute_coverage/CAC_10_007DE.tRNA.position"

def accept_func(feature_name:str):
  return(feature_name.startswith("tRNA:"))
# assert(not accept_func("miRNA:hsa-let-7a-5p"))
# assert(accept_func("tRNA:hsa-let-7a-5p"))

def get_names(feature_name:str):
  full_name = feature_name.split(':')[1]
  parts = full_name.split('-')
  nparts = len(parts)
  transcript = "-".join(parts[0:(nparts-1)])
  anticodon = "-".join(parts[0:(nparts-2)])
  isotype =  "-".join(parts[0:(nparts-3)])
  return(isotype, anticodon, transcript)
# isotype, anticodon, transcript=get_names("tRNA:tRNA-Und-NNN-4-1")
# assert("tRNA-Und" == isotype)
# assert("tRNA-Und-NNN" == anticodon)
# assert("tRNA-Und-NNN-4" == transcript)

def extract_tRNA_position(logger, xml_file, position_file):
  logger.info(f"Reading {xml_file} ...")

  tree = ET.parse(xml_file)
  root = tree.getroot()
  features_nodes = root.find('subjectResult')

  logger.info(f"Parsing {xml_file} ...")
  transcript_map = {}
  query_map = {}
  for feature_group_node in features_nodes.findall('subjectGroup'):
    feature_nodes = feature_group_node.findall('subject')
    for feature_node in feature_nodes:
      feature_name = feature_node.get("name")
      if not accept_func(feature_name):
        continue
      
      isotype, anticodon, transcript = get_names(feature_name)
      if not transcript in transcript_map:
        transcript_map[transcript] = {
          "size":0,
          "isotype":isotype,
          "anticodon":anticodon,
          "queries":{}
        }
      
      q_map = transcript_map[transcript]["queries"]
      
      region_nodes = feature_node.findall('region')
      for region_node in region_nodes:
        t_size = int(region_node.get("size"))      
        transcript_map[transcript]["size"] = max(transcript_map[transcript]["size"], t_size)

        query_nodes = region_node.findall('query')
        for query_node in query_nodes:
          query_name = query_node.get('qname')
          if query_name in q_map:
            continue

          if not query_name in query_map:
            query_map[query_name] = {
              "isotype":{isotype:1},
              "anticodon":{anticodon:1},
              "transcript":{transcript:1},
            }
          else:
            query_map[query_name]["isotype"][isotype] = 1
            query_map[query_name]["anticodon"][anticodon] = 1
            query_map[query_name]["transcript"][transcript] = 1
          
          offset = int(query_node.get('offset'))
          query_count = int(query_node.get('query_count'))
          seq_len = int(query_node.get('seq_len'))

          q_map[query_name] = {
            "offset":offset,
            "query_count":query_count,
            "seq_len":seq_len
          }

  for query_name in query_map.keys():
    q_map = query_map[query_name]
    if len(q_map["transcript"]) == 1:
      q_map["category"] = "Transcript Specific"
    elif len(q_map["anticodon"]) == 1:
      q_map["category"] = "Isodecoder Specific"
    elif len(q_map["isotype"]) == 1:
      q_map["category"] = "Isotype Specific"
    else:
      q_map["category"] = "Not Amino Specific"

  #print(query_map['A00252:306:HLH3MDSX5:1:1127:20238:7827:CLIP_'])
  # for qname, q_map in query_map.items():
  #   print(qname + ":" + q_map['category'])

  logger.info("Calculating position count ...")
  position_map = {}
  for transcript in transcript_map.keys():
    t_size = transcript_map[transcript]["size"]

    position_map[transcript]={
      "Transcript Specific":{
        0:0, 
        t_size-1:0
      },
      "Isodecoder Specific":{
        0:0, 
        t_size-1:0
      },
      "Isotype Specific":{
        0:0, 
        t_size-1:0
      },
      "Not Amino Specific":{
        0:0, 
        t_size-1:0
      }
    }

    q_map = transcript_map[transcript]['queries']
    for qname, v_map in q_map.items():
      category = query_map[qname]["category"]
      offset = v_map["offset"]
      query_count = v_map["query_count"]
      seq_len = v_map["seq_len"]

      c_map = position_map[transcript][category]
      for pos in range(offset, offset+seq_len):
        c_map[pos] = c_map.get(pos, 0) + query_count

  #print(position_map)
  logger.info(f"Saving to file {position_file} ...")
  with open(position_file, "wt") as fout:
    fout.write("Transcript\tCategory\tPosition\tCount\n")
    for transcript in sorted(position_map.keys()):
      t_map = position_map[transcript]
      for cat in sorted(t_map.keys()):
        c_map = t_map[cat]
        for pos in sorted(c_map.keys()):
          fout.write(f"{transcript}\t{cat}\t{pos}\t{c_map[pos]}\n")

  logger.info("Done.")

def main():
  parser = argparse.ArgumentParser(description="Extract tRNA position count.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input xml file')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output position text file")

  args = parser.parse_args()
  
  print(args)
  
  logger = logging.getLogger('tRNA_position')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  extract_tRNA_position(logger, args.input, args.output)

if __name__ == "__main__":
    main()
