rm(list=ls()) 
sample_name='P6121_CP_35'
outFile='P6121_CP_35'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/shengq2/20210616_human_exomeseq/ascat/result/P6121_CP_35')

### Parameter setting end ###

library(ASCAT)
library(data.table)
library(dplyr)

delete_files <- function(pattern){
  temp_files = list.files(pattern = pattern)
  removed=file.remove(temp_files)
}

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
tumourRawBAF_file = paste0(tumourname, '_BAF_rawBAF.txt')

### prepare HTS
print("ascat.prepareHTS...")
ascat.prepareHTS(
  tumourseqfile = tumourseqfile,
  tumourname = tumourname,
  allelecounter_exe = allelecounter_exe,
  alleles.prefix = alleles.prefix,
  loci.prefix = loci.prefix,
  BED_file = BED_file,
  gender = gender, 
  genomeVersion = genomeVersion, 
  nthreads = 8,
  tumourLogR_file = tumourLogR_file,
  tumourBAF_file = tumourBAF_file)
 
delete_files("*.alleleFrequencies*")
system(paste("gzip", shQuote(tumourRawBAF_file)))

### run ASCAT
# gamma=1 (suggested)
print("ascat.loadData...")
ascat.bc = ascat.loadData(Tumor_LogR_file = tumourLogR_file, 
                          Tumor_BAF_file = tumourBAF_file, 
                          gender = gender, 
                          genomeVersion = genomeVersion) 

print("gzip LogR and BAF files ...")
system(paste("gzip", shQuote(tumourLogR_file)))
system(paste("gzip", shQuote(tumourBAF_file)))

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
file.rename(paste0("tumorSep", tumourname, ".png"), paste0(tumourname, ".tumourSep.png"))

print("ascat.aspcf...")
if("penalty" %in% names(my_options)){
  penalty=as.numeric(my_options$penalty)
  cat("penalty=", penalty, "\n")
  ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=gg, penalty=penalty)
}else{
  ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=gg)
}

print("ascat.plotSegmentedData...")
ascat.plotSegmentedData(ascat.bc, img.prefix = paste0(tumourname, ".Segmented."))
file.rename(paste0(tumourname, ".Segmented.", tumourname, ".ASPCF.png"), paste0(tumourname, ".Segmented.ASPCF.png"))

delete_files("*segments*")

print("ascat.runAscat...")
if("min_ploidy" %in% names(my_options) && "max_ploidy" %in% names(my_options)){
  min_ploidy=as.numeric(my_options$min_ploidy)
  cat("min_ploidy=", min_ploidy, "\n")
  max_ploidy=as.numeric(my_options$max_ploidy)
  cat("max_ploidy=", max_ploidy, "\n")
  ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments=TRUE, min_ploidy=min_ploidy, max_ploidy=max_ploidy)
}else{
  ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments=TRUE)
}

delete_files("*.PCFed.txt")

print("ascat.metrics...")
QC = ascat.metrics(ascat.bc, ascat.output)

result=list(ascat.bc=ascat.bc, ascat.output=ascat.output, QC=QC)
saveRDS(result, paste0(tumourname, '.ascat.rds'))

unlink(".cache", recursive = TRUE, force = TRUE)

