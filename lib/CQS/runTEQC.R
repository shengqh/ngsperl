#perform TEQC for bam file target coverage

TEQCprojectName=paste0(outFile,".TEQC")
allBamFiles=read.delim(parSampleFile1,header=FALSE,as.is=TRUE)
targetsfile <- parFile1
#genome="hg38"
library(TEQC)


#target region
print(paste0("Reading target file: ",targetsfile))
targets=get.targets(targetsfile, skip = 0)

#QC for each sample
for (i in 1:nrow(allBamFiles)) {
  bamOne=allBamFiles[i,1]
  sampleName=allBamFiles[i,2]
  if (file.exists(paste0(sampleName,"/index.html"))) {
    print(paste0(sampleName," TEQC result existed. Skip;"))
    next;
  } else {
    print(paste0("Running TEQC on ",sampleName,";"))
  }

  reads= get.reads(bamOne, filetype = "bam")
  if (genome %in% c("hg19","hg38","hg18")) {
    TEQCreport(sampleName = sampleName,targetsName = basename(targetsfile),referenceName = genome, destDir = sampleName,reads =reads,targets = targets,genome = genome)
  } else if (genome=="mm10") {
    TEQCreport(sampleName = sampleName,targetsName = basename(targetsfile),referenceName = genome, destDir = sampleName,reads =reads,targets = targets,genome = "NA",genomesize = 2730871774)
  } else {
    stop(paste0("genome (",genome,") is not supported. Only support hg18, hg19, hg38, mm10"))
  }
}

validResults=allBamFiles[which(file.exists(paste0(allBamFiles[,2],"/index.html"))),2]
print("Running TEQC results merged report;")
multiTEQCreport(singleReportDirs = validResults,
                samplenames = validResults,
                projectName = TEQCprojectName,
                targetsName = basename(targetsfile),
                referenceName = genome,
                destDir = TEQCprojectName)


