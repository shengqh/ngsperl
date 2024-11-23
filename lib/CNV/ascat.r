rm(list=ls()) 
sample_name='P1809_AC_003'
outFile='P1809_AC_003'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/shengq2/20210616_human_exomeseq/ascat/result/P1809_AC_003')

### Parameter setting end ###

library(ASCAT)
library(data.table)
library(dplyr)

opt_tbl=fread(parSampleFile1, header=F) 
tumor_file_map=split(opt_tbl$V1, opt_tbl$V2)
print(tumor_file_map)

tumourname=sample_name
tumourseqfile=tumor_file_map[[tumourname]]

if(!file.exists(tumourseqfile)){
  stop(paste0("File not found: ", tumourseqfile))
}

opt_tbl=fread(parSampleFile2, header=F)
my_options=split(opt_tbl$V1, opt_tbl$V2)
print(my_options)

alleles.prefix=my_options$alleles.prefix
loci.prefix=my_options$loci.prefix
BED_file=my_options$BED_file
GCcontentfile=my_options$GCcontentfile
replictimingfile=my_options$replictimingfile
genomeVersion=my_options$genomeVersion
allelecounter_exe=my_options$allelecounter_exe

opt_tbl=fread(parSampleFile3, header=F)
gender_map=split(opt_tbl$V1, opt_tbl$V2)
gender=gender_map[[tumourname]]

print(paste0("gender=",gender))

tumourLogR_file = paste0(tumourname, '_LogR.txt')
tumourBAF_file = paste0(tumourname, '_BAF.txt')

### prepare HTS
print("ascat.prepareHTS...")
ascat.prepareHTS(
  tumourseqfile = tumourseqfile,
  tumourname = tumourname,
  allelecounter_exe = 'alleleCounter',
  alleles.prefix = alleles.prefix,
  loci.prefix = loci.prefix,
  BED_file = BED_file,
  gender = gender, 
  genomeVersion = genomeVersion, 
  nthreads = 8,
  tumourLogR_file = tumourLogR_file,
  tumourBAF_file = tumourBAF_file)
 
### run ASCAT
# gamma=1 (suggested)
print("ascat.loadData...")
ascat.bc = ascat.loadData(Tumor_LogR_file = tumourLogR_file, 
                          Tumor_BAF_file = tumourBAF_file, 
                          gender = gender, 
                          genomeVersion = genomeVersion) 

print("ascat.plotRawData before correlation...")                          
ascat.plotRawData(ascat.bc, img.prefix = paste0(tumourname, ".Before_correction."))
file.rename(paste0(tumourname, ".Before_correction.", tumourname, ".tumour.png"), paste0(tumourname, ".Before_correction.tumour.png"))

print("ascat.correctLogR...")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = GCcontentfile, replictimingfile = replictimingfile)

print("ascat.plotRawData after correlation...")
ascat.plotRawData(ascat.bc, img.prefix = paste0(tumourname, ".After_correction."))
file.rename(paste0(tumourname, ".After_correction.", tumourname, ".tumour.png"), paste0(tumourname, ".After_correction.tumour.png"))

print("ascat.predictGermlineGenotypes...")
gg = ascat.predictGermlineGenotypes(ascat.bc, platform = "WGS_hg38_50X")

print("ascat.aspcf...")
ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=gg)

print("ascat.plotSegmentedData...")
ascat.plotSegmentedData(ascat.bc, img.prefix = paste0(tumourname, ".Segmented."))
file.rename(paste0(tumourname, ".Segmented.", tumourname, ".ASPCF.png"), paste0(tumourname, ".Segmented.ASPCF.png"))

print("ascat.runAscat...")
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments=TRUE)

print("ascat.metrics...")
QC = ascat.metrics(ascat.bc, ascat.output)

result=list(ascat.bc=ascat.bc, ascat.output=ascat.output, QC=QC)
saveRDS(result, paste0(tumourname, '.ascat.rds'))

#remove temporary files
for(pat in c("*alleleFrequencies*", "*.PCFed.txt", "*_rawBAF*", "*segments*")){
  temp_files = list.files(pattern = pat)
  removed=file.remove(temp_files)
}

file.rename(paste0("tumorSep", tumourname, ".png"), paste0(tumourname, ".tumourSep.png"))

system(paste("gzip", shQuote(tumourLogR_file)))
system(paste("gzip", shQuote(tumourBAF_file)))

unlink(".cache", recursive = TRUE, force = TRUE)
