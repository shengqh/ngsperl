import subprocess
import os.path
import re
import argparse

parser = argparse.ArgumentParser(description="annovate splicing with protein position by annovar.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input annovar _multianno.txt file', required=True)
parser.add_argument('-d', '--database_folder', action='store', nargs='?', help='Annovar database folder', required=True)
parser.add_argument('-b', '--buildver', action='store', nargs='?', help='Annovar buildver', required=True)
parser.add_argument('-s', '--splicing_threshold', action='store', nargs='?', default=2, help='Annovar splicing threshold (default=2)')
parser.add_argument('-o', '--output', action='store', nargs='?', help='Output annovar _multianno.txt file', required=True)

args = parser.parse_args()

inputfile=args.input
annovar_db=args.database_folder
outputfile=args.output
annovar_buildver=args.buildver
splicing_distance=int(args.splicing_threshold)

# inputfile="/scratch/cqs/shengq1/dnaseq/20160829_liuqi_gene_panel/bwa_refine_hc_gvcf_vqsr_annovar/result/liuqi_gene/liuqi_gene.pass.annovar.hg19_multianno.txt"
# annovar_db="/scratch/cqs/shengq1/references/annovar/humandb/"
# outputfile="/scratch/cqs/shengq1/dnaseq/20160829_liuqi_gene_panel/bwa_refine_hc_gvcf_vqsr_annovar/result/liuqi_gene/liuqi_gene.final.annovar.hg19_multianno.txt"
# annovar_buildver="hg19"
# splicing_distance=2

possible_positions = [j for j in range(-splicing_distance, splicing_distance + 1) if j != 0]

with open(inputfile, 'r') as f:
  headers = f.readline().rstrip().split('\t')
  funcRefGeneIndex=headers.index("Func.refGene")
  geneDetailRefGeneIndex = headers.index("GeneDetail.refGene")
  aachangeRefGeneIndex=headers.index("AAChange.refGene")

annovarinputfile = outputfile + ".avinput"
with open(inputfile, 'r') as f:
  with open(annovarinputfile, 'w') as snvw:
    for line in f:
      parts = line.rstrip('\r\n').split('\t')
      if parts[funcRefGeneIndex] != "splicing":
        continue
      
      chro = parts[0]
      start = int(parts[1])
      for pos in possible_positions:
        snvw.write("%s\t%d\t%d\tA\tG\t%d\n" % (chro, start + pos, start + pos, start))

args=["table_annovar.pl",
      annovarinputfile,
      annovar_db,
      "-buildver",
      annovar_buildver,
      "-protocol",
      "refGene",
      "-operation",
      "g",
      "--outfile",
      outputfile,
      "--remove",
      "--otherinfo"]
     
subprocess.call(args)

annovar_outputfile = outputfile + "." + annovar_buildver + "_multianno.txt"

if os.path.isfile(annovar_outputfile):
  splicing_map = {}
  prog = re.compile("p\.\w(\d+)[\w|\?]")
  with open(annovar_outputfile, "r") as f:
    splicingHeaders = f.readline().rstrip().split('\t')
    splicingFuncRefGeneIndex=splicingHeaders.index("Func.refGene")
    splicingAAChangeRefGeneIndex=splicingHeaders.index("AAChange.refGene")
    
    for line in f:
      parts = line.rstrip('\r\n').split('\t')
      funcRefGene = parts[splicingFuncRefGeneIndex]
      if(funcRefGene == "splicing" or funcRefGene == "intronic"):
        continue
      
      chrom = parts[0]
      originalposition = int(parts[-1])
      position = int(parts[1])
      distance = abs(position - originalposition)
      if funcRefGene != "exonic":
        anno = funcRefGene
      else:
        anno = {}
        for a in parts[splicingAAChangeRefGeneIndex].split(','):
          if a != 'UNKNOWN':
            aparts = a.split(':')
            ma = prog.match(aparts[-1])
            if ma is None:
              print(line)
            else:
              anno[aparts[1]] = ':'.join(aparts[0:3]) + ':p.X' + ma.group(1) + 'X'
      
      locus = parts[0] + ":" + str(originalposition)
      if locus in splicing_map:
        oldone = splicing_map[locus]
        if oldone[0] > distance:
          splicing_map[locus] = [distance, funcRefGene, anno]
      else:
        splicing_map[locus] = [distance, funcRefGene, anno]

  outputTemp = outputfile + ".tmp"
  with open(inputfile, 'r') as f:
    with open(outputTemp, 'w') as w:
      for line in f:
        parts = line.rstrip('\r\n').split('\t')
        if parts[funcRefGeneIndex] != "splicing":
          w.write(line)
          continue
        
        locus=  parts[0] + ":" + parts[1]
        if locus not in splicing_map:
          w.write(line)
          continue
        
        values = splicing_map[locus]
        if values[1] != 'exonic':
          parts[aachangeRefGeneIndex] = values[2]
        else:
          anno = values[2]
          lst = []
          for detail in parts[geneDetailRefGeneIndex].split(','):
            detailparts = detail.split(':')
            trans = detailparts[0]
            if trans in anno :
              lst.append(anno[trans])
            else:
              lst.append('')
          parts[aachangeRefGeneIndex] = ','.join(lst)
          
        w.write("%s\n" % ("\t".join(parts)))
  os.remove(annovarinputfile)
  os.remove(annovar_outputfile)

  if os.path.isfile(outputfile):
    os.remove(outputfile)
  os.rename(outputTemp, outputfile)
  
print("annotate splicing by annovar done.")
