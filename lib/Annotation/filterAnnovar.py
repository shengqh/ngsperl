import argparse
import re

parser = argparse.ArgumentParser(description="filter Annovar result for truncating/nonsense SNVs with/without ExAC/1000G limitation",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

DEBUG = False
NOT_DEBUG = not DEBUG

parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help='Input annovar result file from NGSPERL')
parser.add_argument('-e', '--threshold', action='store', nargs='?', default=1.0, help='Maximum ExAC/1000G allele frequency (default=1.0, no filter)')
parser.add_argument('-r', '--sample_name_pattern', action='store', nargs='?', default="", help='Sample name regex pattern for extraction those samples only')
parser.add_argument('-o', '--output_prefix', action='store', nargs='?', required=NOT_DEBUG, help='Output file prefix')

args = parser.parse_args()

print(args)

inputfile = args.input
outputprefix = args.output_prefix
maxExacValue = float(args.threshold)
sampleNamePattern = args.sample_name_pattern

if DEBUG:
  inputfile = "H:/shengquanhu/projects/20161013_liuqi_gene_panel/bwa_refine_hc_gvcf_vqsr_annovar/Adenoma.pass.annovar.final.tsv"
  outputprefix = "H:/shengquanhu/projects/20161013_liuqi_gene_panel/bwa_refine_hc_gvcf_vqsr_annovar/Adenoma.pass.annovar.final.exac0.01"
  maxExacValue = 0.01

accepted = ["frameshift deletion", "frameshift insertion", "frameshift substitution", "nonsynonymous SNV", "stopgain", "stoploss"]

def getKey(item): return item[0]
def getDicValueCount(item): return len(item[1])


filtered = []
genes = {}
with open(outputprefix + ".filtered.tsv", 'w') as sw:
  with open(inputfile, 'r') as f:
    for line in f:
      if(not line.startswith("#")):
        break
      sw.write(line)

    headers = line.rstrip().split()
    formatIndex = headers.index("FORMAT")

    snvHeaderIndecies = range(0, formatIndex)
    sampleIndecies = range(formatIndex + 1, len(headers))
    if sampleNamePattern != "":
      sampleIndecies = [index for index in sampleIndecies if re.search(sampleNamePattern, headers[index])]

    sw.write("%s\t%s\t%s\n" %("\t".join(headers[i] for i in snvHeaderIndecies), headers[formatIndex], "\t".join(headers[i] for i in sampleIndecies)))

    funcIndex = headers.index("Func.refGene")
    geneIndex = headers.index("Gene.refGene")
    refIndex = headers.index("ExonicFunc.refGene")
    exacIndex = headers.index("ExAC_ALL")
    g1000Index = headers.index("1000g2015aug_all")
    sampleCount = len(sampleIndecies)

    print("sampleCount=%d" % (sampleCount))
    if sampleCount == 0:
      raise ValueError("Sample count is zero, check your sample name regex pattern: " + sampleNamePattern)

    for line in f:
      parts = line.rstrip().split('\t')
      gene = parts[geneIndex]
      if(parts[funcIndex] == "splicing" or parts[refIndex] in accepted):
        if (exacIndex != -1) and (parts[exacIndex] != "") and float(parts[exacIndex]) > maxExacValue:
          continue
          
        if (g1000Index != -1) and (parts[g1000Index] != "") and float(parts[g1000Index]) > maxExacValue:
          continue
          
        freq = len([idx for idx in sampleIndecies if parts[idx].startswith("0/1") or parts[idx].startswith("1/1")])
        if freq == 0:
          continue

        sw.write("%s\t%s\t%s\n" %("\t".join(parts[i] for i in snvHeaderIndecies), parts[formatIndex], "\t".join(parts[i] for i in sampleIndecies )))

        for idx in sampleIndecies:
          if(parts[idx].startswith("0/1") or parts[idx].startswith("1/1")):
            parts[idx] = "1"
          else:
            parts[idx] = "0"

        freqperc = float(freq) / sampleCount
        freqfold = "NA" if exacIndex == -1 or not parts[exacIndex] else "100" if float(parts[exacIndex]) == 0 else str(freqperc / float(parts[exacIndex]))
        filtered.append([freqperc, "%s\t%f\t%s\t\t%s\n" % ("\t".join(parts[i] for i in snvHeaderIndecies), freqperc, freqfold, "\t".join(parts[i] for i in sampleIndecies))])

        curgene = parts[geneIndex]
        if len(curgene) > 0:
          if curgene in genes:
            curgenemap = genes[curgene]
            for idx in sampleIndecies:
              if parts[idx] == "1":
                curgenemap[idx] = "1"
          else:
            curgenemap = {}
            genes[curgene] = curgenemap
            for idx in sampleIndecies:
              curgenemap[idx] = parts[idx]

    fsorted = sorted(filtered, key=getKey, reverse=True)
    with open(outputprefix + ".snv.tsv", 'w') as snvw:
      snvw.write("%s\tFrequency\tFrequencyFoldChange\tFormat\t%s\n" % ("\t".join(headers[i] for i in snvHeaderIndecies), "\t".join(headers[i] for i in sampleIndecies)))
      for d in fsorted:
        snvw.write(d[1])

    genelist = []
    for gene in genes:
      genemap=genes[gene]
      freq = [v for k,v in genemap.items() if v == "1"]
      freqperc = len(freq) * 1.0 / sampleCount
      genelist.append([freqperc, gene])
    gsorted = sorted(genelist, key=getKey, reverse=True)
    with open(outputprefix + ".gene.tsv", 'w') as genew:
      genew.write("Gene\tFrequency\t%s\n" % ("\t".join(headers[i] for i in sampleIndecies)))
      for d in gsorted:
        genemap=genes[d[1]]
        genew.write("%s\t%f\t%s\n" %(d[1], d[0], "\t".join(genemap[i] for i in sampleIndecies )))

print("done.")
