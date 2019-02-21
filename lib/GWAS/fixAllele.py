import argparse
import sys
import logging
import os
from Bio import SeqIO
import subprocess

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Fix allele in plink bim file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input plink file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output plink file", required=NotDEBUG)
parser.add_argument('-b', '--manifest', action='store', nargs='?', help="Manifest file in csv format", required=NotDEBUG)
parser.add_argument('-f', '--fasta', action='store', nargs='?', help="Genome fasta file", required=NotDEBUG)

args = parser.parse_args()

if DEBUG:
  args.input = "H:/shengquanhu/projects/macrae_linton/2118-JB_GSProject"
  args.output = "H:/shengquanhu/projects/macrae_linton/2118-JB_GSProject.fixed"
  args.manifest = "H:/shengquanhu/projects/macrae_linton/MEGAEx_BioVU_15075710_A1.csv"
  args.fasta = "H:/shengquanhu/projects/database/hg19/hg19_16569_M.fa"

logger = logging.getLogger('fixAllele')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

baseComps = {'A':'T',
             'T':'A',
             'G':'C',
             'C':'G',
             'N':'N',
             '0':'0'}

name_index = 0
position_index = 1
illumina_allele1_index = 2
illumina_allele2_index = 3
alt_allele_index = 4
ref_allele_index = 5

def readSnpFile(refFile):
  snpMap = {}
  with open(refFile, "r") as fin:
    fin.readline()
    for line in fin:
      parts = line.rstrip().split('\t')
      chr = parts[0]
      position = int(parts[1])
      name = parts[2]
      if chr not in snpMap:
        chrSnpMap = {}
        snpMap[chr] = chrSnpMap
      else:
        chrSnpMap = snpMap[chr]
      chrSnpMap[name] = [name, position, parts[5], parts[6], parts[4], parts[3]]
  return(snpMap)
        
def writeSnpFile(refFile, snpMap):
  with open(refFile, "w") as fout:
    fout.write("chr\tposition\tsnp\tref_allele\talt_allele\tillumina_allele1\tillumina_allele2\n")
    for chr in sorted(snpMap):
      chrSnpMap = snpMap[chr]
      snpInfos = [snpInfo for snp, snpInfo in chrSnpMap.iteritems()]
      snpInfos.sort(key=lambda s:s[position_index])
      for snpInfo in snpInfos:
        fout.write(chr + 
                   "\t" + str(snpInfo[position_index]) + 
                   "\t" + snpInfo[name_index] + 
                   "\t" + snpInfo[ref_allele_index] +
                   "\t" + snpInfo[alt_allele_index] + 
                   "\t" + snpInfo[illumina_allele1_index] + 
                   "\t" + snpInfo[illumina_allele2_index] + "\n")
    
refFile = args.manifest + ".ref.txt"
if os.path.isfile(refFile):
  logger.info("Reading " + refFile + "...")
  snpMap = readSnpFile(refFile)
else:
  logger.info("Reading " + args.manifest + "...")
  snpMap = {}
  with open(args.manifest, "r") as fh:
    for line in fh:
      if line.startswith("[Assay]"):
        headerline = fh.next()
        headers = headerline.split(",")
        nameIndex = headers.index("Name")
        snpIndex = headers.index("SNP")
        chrIndex = headers.index("Chr")
        positionIndex = headers.index("MapInfo")
        break
    
    for line in fh:
      parts = line.split(',')
      name = parts[nameIndex]
      snp = parts[snpIndex]
      allele1 = snp[1]
      allele2 = snp[3]
      chr = parts[chrIndex]
      position = int(parts[positionIndex])
      if chr not in snpMap:
        snpMap[chr]  = {}
      snpMap[chr][name] = [name, position, allele1, allele2, '', ''] #name, position, illumina_allele1, illumina_allele2, ref_allele, alt_allele
  
  logger.info("Reading " + args.fasta + "...")
  with open(args.fasta, "r") as fin:    
    for record in SeqIO.parse(fin,'fasta'):
      refChr = record.id
      logger.info("Processing " + refChr + " ...")
      
      if refChr.startswith("chr"):
        refChr = refChr[3:]
        
      if (refChr not in snpMap) and (refChr == "M"):
        refChr = "MT"
      
      if record.id in snpMap:
        chrSnpMap = snpMap[record.id]
        for snp, snpInfo in chrSnpMap.iteritems():
          if snpInfo[illumina_allele1_index] == 'I' or snpInfo[illumina_allele1_index] == 'D':
            continue
          
          refBase = record.seq[snpInfo[position_index]].upper()
          snpInfo[ref_allele_index] = refBase
          if snpInfo[illumina_allele1_index] == refBase:
            snpInfo[alt_allele_index] = snpInfo[illumina_allele2_index]
            continue
          
          if snpInfo[illumina_allele2_index] == refBase:
            snpInfo[alt_allele_index] = snpInfo[illumina_allele1_index]
            continue
          
          allele1comp = baseComps[snpInfo[illumina_allele1_index]]
          allele2comp = baseComps[snpInfo[illumina_allele2_index]]
          if allele1comp == refBase:
            snpInfo[alt_allele_index] = allele2comp
            continue
          
          if allele2comp == refBase:
            snpInfo[alt_allele_index] = allele1comp
            continue
  
  logger.info("Writing " + refFile + "...")
  writeSnpFile(refFile, snpMap)

expectBimFile = args.output + ".expect.bim"
with open(args.input + ".bim", "r") as fin:
  with open(expectBimFile, "w") as fout_bim:
    for line in fin:
      parts = line.rstrip().split('\t')
      chr = parts[0]
      name = parts[1]
      position = int(parts[3])
      minor_allele = parts[4]
      major_allele = parts[5]
       
      if chr not in snpMap:
        continue
   
      chrSnpMap = snpMap[chr]
      if name not in chrSnpMap:
        continue
  
      snpInfo = chrSnpMap[name]
      if snpInfo[alt_allele_index] == '':
        continue
      
      vdata = "%s\t%s\t0\t%d\t%s\t%s\t%s\t%s" % (chr, name, position, snpInfo[alt_allele_index], snpInfo[ref_allele_index], minor_allele, major_allele)
      if major_allele == '0':
        fout_bim.write(vdata + "\tFill\n")
        continue
      
      if major_allele == snpInfo[ref_allele_index]:
        if minor_allele == '0':
          fout_bim.write(vdata + "\tFill\n")
        elif minor_allele != snpInfo[alt_allele_index]:
          logger.error("Conflict: " + vdata)
        continue

      if major_allele == snpInfo[alt_allele_index]:
        if (minor_allele == snpInfo[ref_allele_index] or minor_allele == '0'):
          fout_bim.write(vdata + "\tSwitch\n")
        else:
          logger.error("Conflict: " + vdata)
        continue
      
      major_comp = baseComps[major_allele]
      minor_comp = baseComps[minor_allele]
      if major_comp == snpInfo[ref_allele_index]:
        if (minor_comp == snpInfo[alt_allele_index] or minor_allele == '0'):
          fout_bim.write(vdata + "\tFlip\n")
        else:
          logger.error("Conflict: " + vdata)
        continue
      
      if major_comp == snpInfo[alt_allele_index]:
        if (minor_comp == snpInfo[ref_allele_index] or minor_allele == '0'):
          fout_bim.write(vdata + "\tFlip+Switch\n")
        else:
          logger.error("Conflict: " + vdata)
        continue

flipListFile = args.output + ".flip.list"
flipFile = args.output + ".flip"
refAlleleFile = args.output + ".flip.refAllele"


os.system("grep Flip " + expectBimFile + " | cut -f2 > " + flipListFile)
os.system("plink --bfile " + args.input + " --flip " + flipListFile + " --make-bed --out " + flipFile)
os.system("plink --bfile  " + flipFile + " --keep-allele-order --real-ref-alleles --a2-allele " + expectBimFile + " 6 2 '#' --make-bed --out " + refAlleleFile)
os.system("plink --bfile  " + refAlleleFile + " --keep-allele-order --real-ref-alleles --a1-allele " + expectBimFile + " 5 2 '#' --make-bed --out " + args.output)

if os.path.isfile(args.output + ".bed"):
  os.system("rm " + flipFile + ".bed " + flipFile + ".bim " + flipFile + ".fam " + flipFile + ".hh " + refAlleleFile + ".bed " + refAlleleFile + ".bim " + refAlleleFile + ".fam " + refAlleleFile + ".hh ")

logger.info("Done.")
