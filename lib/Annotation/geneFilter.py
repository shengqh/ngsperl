import argparse
import re

parser = argparse.ArgumentParser(description="Export filtered genotype information of gene", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

DEBUG = False
NOT_DEBUG = not DEBUG

parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help='Input filtered annovar result file from FilterAnnovar')
parser.add_argument('-g', '--genes', action='store', nargs='?', required=NOT_DEBUG, help='Gene name list, seperated by ","')
parser.add_argument('-o', '--output', action='store', nargs='?', required=NOT_DEBUG, help='Output file')

args = parser.parse_args()

print(args)

inputfile = args.input
genes = args.genes
outputfile = args.output

if DEBUG:
  inputfile = "H:/temp/Adenoma.base1.exac0.001.filtered.tsv"
  outputfile = "H:/temp/Adenoma.base1.exac0.001.filtered.gene.snv.tsv"
  genes = "LRDR"

genelist = genes.split(",")

with open(outputfile, 'w') as sw:
  bFirst = True
  for gene in genelist:
    with open(inputfile, 'r') as f:
      for line in f:
        if(not line.startswith("#")):
          break

      if bFirst:
        headers = line.rstrip().split()
        formatIndex = headers.index("FORMAT")
        geneIndex = headers.index("Gene.refGene")
        aaIndex = headers.index("AAChange.refGene")
        snvHeaderIndecies = range(0, formatIndex)
        sampleIndecies = range(formatIndex + 1, len(headers))
        sampleCount = len(sampleIndecies)
        print("sampleCount=%d" % (sampleCount))
        if sampleCount == 0:
          raise ValueError("Sample count is zero, check your file : " + inputfile)
        
        sw.write(line)
        bFirst = False

      for line in f:
        parts = line.rstrip().split('\t')
        curGenes = parts[geneIndex].split(",")
        if gene in curGenes:
          parts[aaIndex] = parts[aaIndex].split(",")[0]
          sw.write("%s" % "\t".join([parts[idx] for idx in snvHeaderIndecies]))
          for idx in sampleIndecies:
            if parts[idx].startswith("0/0"):
              sw.write("\tAA")
            elif parts[idx].startswith("0/1"):
              sw.write("\tAB")
            elif parts[idx].startswith("1/1"):
              sw.write("\tBB")
            else:
              sw.write("\t")
          sw.write("\n")

print("done.")
