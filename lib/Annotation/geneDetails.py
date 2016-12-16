import argparse
import re

parser = argparse.ArgumentParser(description="Export detail genotype information of gene", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

DEBUG = True
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
  outputfile = "H:/temp/Adenoma.base1.exac0.001.filtered.gene.details.tsv"
  genes = "APC"

genelist = genes.split(",")
with open(outputfile, 'w') as sw:
  sw.write("Gene_Symbol\tt_alt_count\tt_ref_count\tTot_Count\tAlt%\tChromosome\tStart_position\tVariant_Classification\tProtein_Change\tTumor_Sample_UUID\n")
  for gene in genelist:
    with open(inputfile, 'r') as f:
      for line in f:
        if(not line.startswith("#")):
          break

      headers = line.rstrip().split()
      formatIndex = headers.index("FORMAT")

      geneIndex = headers.index("Gene.refGene")
      funcIndex = headers.index("Func.refGene")
      refIndex = headers.index("ExonicFunc.refGene")
      aaIndex = headers.index("AAChange.refGene")

      snvHeaderIndecies = range(0, formatIndex)
      sampleIndecies = range(formatIndex + 1, len(headers))

      sampleCount = len(sampleIndecies)
      print("sampleCount=%d" % (sampleCount))
      if sampleCount == 0:
        raise ValueError("Sample count is zero, check your file : " + inputfile)

      for line in f:
        parts = line.rstrip().split('\t')
        curGenes = parts[geneIndex].split(",")
        if gene in curGenes:
          for idx in sampleIndecies:
            print("%d" % idx)
            if parts[idx].startswith("0/1") or parts[idx].startswith("1/1"):
              counts = parts[idx].split(":")[1].split(",")
              refCount=int(counts[0])
              altCount=int(counts[1])
              totalCount=refCount+altCount
              frequency=altCount*1.0/totalCount
              variantClass = parts[refIndex] if parts[refIndex] != "" else parts[funcIndex]
              aachange = parts[aaIndex]
              if aachange != "":
                aachange = ",".join(aa for aa in aachange.split(",") if aa.split(":")[0] == gene)
              sw.write("%s\t%d\t%d\t%d\t%f\t%s\t%s\t%s\t%s\t%s\n" % (gene, altCount, refCount, totalCount, frequency, parts[0], parts[1], variantClass, aachange, headers[idx] ))

print("done.")
