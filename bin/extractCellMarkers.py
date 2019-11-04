import argparse
import sys
import logging
import os
import csv
import re

DEBUG=True
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Extract cell markers from CellMarker database",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input CellMarker file', required=NotDEBUG)
parser.add_argument('-t', '--ontology', action='store', nargs='?', help='Input cell ontology file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output Output file", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  #args.input = "H:/shengquanhu/projects/database/scRNA/Mouse_cell_markers.txt"
  args.ontology = "H:/shengquanhu/projects/database/scRNA/cl-base.obo"
  #args.output = "H:/shengquanhu/projects/database/scRNA/Mouse_cell_markers.filered.txt"
  args.input = "Z:/shengq1/scRNA/Human_cell_markers.txt"
  args.output = "Z:/shengq1/scRNA/Human_cell_markers.filered.txt"

logger = logging.getLogger('cellMarkers')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')


def simplify(marker_genes):
    result = ""
    while '[' in marker_genes:
        parts = marker_genes.split('[', 1)
        result = result + parts[0]
        items = parts[1].split(']', 1)
        inter_genes = items[0]
        marker_genes = items[1]
        inter_genes_parts = inter_genes.split(',')
        result = result + inter_genes_parts[0]
    result = result + marker_genes
    return result

ontologyMap = {}
with open(args.ontology, 'rt') as fin:
  for line in fin:
    if line.startswith("[Term]"):
      clId = ""
      clName = ""
    elif line.startswith("id: CL:"):
      clId = line[7:].rstrip()
    elif line.startswith("name: ") and clId != "":
      clName = line[6:].rstrip()
      ontologyMap["CL_" + clId] = clName
      
ontologyGeneMap = {}
with open(args.input) as tsvfile:
  reader = csv.DictReader(tsvfile, dialect='excel-tab')
  with open(args.output, 'wt') as fout:
    for row in reader:
      olo = row["CellOntologyID"]
      
      if olo == "NA":
        continue
      
      if olo not in ontologyMap:
        print("Cannot find ontology " + olo)
        continue
      
      genes = row["geneSymbol"]
      
      if not olo in ontologyGeneMap:
        ontologyGeneMap[olo] = {}
      
      geneCountMap = ontologyGeneMap[olo]
      
      genes = simplify(genes)
      
      geneparts = genes.split(',')
      for gene in geneparts:
        gene = gene.strip()
        if gene != 'NA':
          if not gene in geneCountMap:
            geneCountMap[gene] = 1
          else:
            geneCountMap[gene] = geneCountMap[gene] + 1

sorted_cell = sorted(ontologyGeneMap.keys(), key=lambda x: ontologyMap[x])
for key in sorted_cell:
  sorted_d = sorted(ontologyGeneMap[key].items(), key=lambda x: x[1], reverse=True)
  print("%s : %s : %s" % (key, ontologyMap[key], sorted_d))

with open(args.output, "wt") as fout:
  fout.write("cellName\tcellId\tmarkerGenes\tmarkerGeneRefCount\n")
  for key in sorted_cell:
    sorted_d = sorted(ontologyGeneMap[key].items(), key=lambda x: x[1], reverse=True)
    fout.write("%s\t%s\t%s\t%s\n" % (ontologyMap[key], key, ','.join([v[0] for v in sorted_d]), ','.join([str(v[1]) for v in sorted_d])))
            
logger.info("done.")

