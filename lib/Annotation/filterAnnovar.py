import argparse
import logging

parser = argparse.ArgumentParser(description="filter Annovar by truncating/nonsense SNV with/without ExAC annotation",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

DEBUG=False
NOT_DEBUG=not DEBUG

parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help='Input annovar result file from NGSPERL')
parser.add_argument('-e', '--exac_threshold', action='store', nargs='?', default=1.0, help='Maximum ExAC value (default=1.0, no filter)')
parser.add_argument('-o', '--output', action='store', nargs='?', required=NOT_DEBUG, help='Output annovar result file')

args = parser.parse_args()

inputfile=args.input
outputfile=args.output
maxExacValue=float(args.exac_threshold)

if DEBUG:
  inputfile="/scratch/cqs/shengq1/dnaseq/20161013_liuqi_gene_panel/bwa_refine_hc_gvcf_vqsr_annovar/result/Adenoma/Adenoma.pass.annovar.final.tsv"
  outputfile="/scratch/cqs/shengq1/dnaseq/20161013_liuqi_gene_panel/bwa_refine_hc_gvcf_vqsr_annovar/result/Adenoma/Adenoma.pass.annovar.final.exac0.01.tsv"
  maxExacValue=1.0

accepted = ["frameshift deletion", "frameshift insertion","frameshift substitution","nonsynonymous SNV","stopgain","stoploss"]

def getKey(item): return item[0]
def getDicValueCount(item): return len(item[1])

with open(inputfile, 'r') as f:
  for line in f:
    if(not line.startswith("#")):
      break

  headers = line.rstrip().split()
  formatIndex = headers.index("FORMAT")

  snvHeaderIndecies = range(0, formatIndex)
  sampleIndecies = range(formatIndex+1,len(headers))

  funcIndex=headers.index("Func.refGene")
  geneIndex = headers.index("Gene.refGene")
  refIndex=headers.index("ExonicFunc.refGene")
  exacIndex=headers.index("ExAC_ALL")
  sampleCount = len(sampleIndecies)
  
  print("sampleCount=%d" % (sampleCount))

  with open(outputfile, 'w') as snvw:
    snvw.write("%s\tFrequency\tFrequencyFoldChange\tFormat\t%s\n" % ("\t".join(headers[i] for i in snvHeaderIndecies), "\t".join(headers[i] for i in sampleIndecies)))
    filtered = []
    for line in f:
      parts = line.split('\t')
      gene = parts[geneIndex]
      if(parts[funcIndex] == "splicing" or parts[refIndex] in accepted):
        if(exacIndex == -1 or not parts[exacIndex] or float(parts[exacIndex]) < maxExacValue):
          freq = 0.0
          for idx in sampleIndecies:
            if(parts[idx].startswith("0/1") or parts[idx].startswith("1/1")):
              parts[idx] = "1"
              freq += 1.0
            else:
              parts[idx] = "0"
          freqperc = freq / sampleCount
          freqfold = "NA" if exacIndex == -1 or not parts[exacIndex] else "100" if float(parts[exacIndex])== 0 else str(freqperc / float(parts[exacIndex]))
          filtered.append([freqperc, "%s\t%f\t%s\t\t%s\n" % ("\t".join(parts[i] for i in snvHeaderIndecies),freqperc, freqfold, "\t".join(parts[i] for i in sampleIndecies))])
    fsorted=sorted(filtered, key=getKey, reverse=True)
    for d in fsorted:
      snvw.write(d[1])
  f.close()
     
print("all done.")