rm(list=ls()) 
outFile='AK6383'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose_edgeR_inCluster_bySample/result/AK6383.edgeR.files.csv'
parFile2=''
parFile3=''


setwd('/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose_edgeR_inCluster_bySample_WebGestalt/result')

### Parameter setting end ###

source("WebGestaltReportFunctions.r")
options(bitmapType='cairo')

library(WebGestaltR)

sigGeneList<-read.csv(parFile1)
sigGeneList$sigGenenameFile<-paste0(dirname(parFile1), "/", sigGeneList$sigGenenameFile)

params_def=read.table(parSampleFile1, sep="\t", stringsAsFactors = F)
params <- setNames(as.character(params_def$V1), params_def$V2)

organism=params['organism']
interestGeneType=params['interestGeneType']
referenceSet=params['referenceSet']
outputDirectory = getwd()

cat("organism=", organism, "\n")
cat("outputDirectory=", outputDirectory, "\n")

if(interestGeneType == ""){
  interestGeneType="genesymbol"
}

if(referenceSet == ""){
  referenceSet="genome"
}
cat("interestGeneType=", interestGeneType, "\n")
cat("referenceSet=", referenceSet, "\n")

allres=NULL
idx=1
for(idx in c(1:nrow(sigGeneList))) {
  sampleName = sigGeneList$comparison[idx]
  geneFile = sigGeneList$sigGenenameFile[idx]
  
  cat("sampleName=", sampleName, "\n")
  cat("geneFile=", geneFile, "\n")

  info = file.info(geneFile)
  if (all(is.na(info$size))) {
    stop(paste0("Gene file is not exist: ", geneFile, "\n"))
  }
  
  if (all(info$size == 0)) {
    error_msg = paste0("Gene file is empty: ", geneFile, "\n")
    cat(error_msg)
    writeLines(error_msg, paste0(outputDirectory, "/", basename(geneFile), ".isempty"))
    next
  }
  
  if (grepl(".csv$",basename(geneFile))) { #csv file
    geneList<-read.csv(geneFile,header=TRUE,stringsAsFactors=FALSE)
  } else {
    geneList<-read.table(geneFile,header=TRUE,sep="\t",stringsAsFactors=FALSE)
  }
  
  if (ncol(geneList)==1) {
    genes<-readLines(geneFile)
  } else { #try to find gene column
    geneCol=getGeneCol(geneList)[["colName"]]
    genes<-geneList[,geneCol]
  }
  genes=unique(genes)

  if(length(genes) < 10){  
    error_msg = paste0("Less than 10 DE genes: ", geneFile, "\n")
    cat(error_msg)
    writeLines(error_msg, paste0(outputDirectory, "/", basename(geneFile), ".less10"))
    next
  }

  referenceGeneFile = gsub('.sig_genename.txt$', '.csv', geneFile)
  if(file.exists(referenceGeneFile)){
    referenceGene=rownames(read.csv(referenceGeneFile, row.names=1))
  }else{
    referenceGene=NULL
  }
  #print(referenceGene)

  gName=gsub(".sig_genename.txt","",basename(geneFile))

  enrichDatabases<-c("geneontology_Biological_Process_noRedundant", 
                     "geneontology_Cellular_Component_noRedundant", 
                     "geneontology_Molecular_Function_noRedundant",
                     "pathway_KEGG" 
                     #"pathway_Wikipathway", 
                     #"network_miRNA_target",
                     #"network_PPI_BIOGRID", 
                     #"network_Transcription_Factor_target"
  )
  names(enrichDatabases)<-c("GO_BP_nr","GO_CC_nr","GO_MF_nr", "KEGG")
  
  for(dbName in names(enrichDatabases)){
    enrichDatabase = enrichDatabases[dbName]
    projectName=paste0(gName, "_", dbName)
    subdir = paste0(outputDirectory, "/Project_", projectName)
    resultFile=paste0(subdir, "/enrichment_results_", projectName, ".txt")
    allres=rbind(allres, data.frame("File"=resultFile, "Database"=enrichDatabase, "Comparison"=sampleName))
    if(file.exists(resultFile)){
      next
    }
    temp=WebGestaltR(enrichMethod="ORA",
                     organism=organism,
                     enrichDatabase=enrichDatabase,
                     interestGene=genes,
                     interestGeneType=interestGeneType,
                     referenceGene=referenceGene,
                     referenceGeneType=interestGeneType,
                     referenceSet=referenceSet,
                     isOutput=TRUE,
                     minNum=5,
                     outputDirectory=outputDirectory,
                     projectName=projectName)
    if (is.null(temp)) { #no enrichment. report top 5 categories
      warning(paste0("No significant category (FDR<=0.05) identified in ",enrichDatabase,", reporting top 10 categories instead."))
      temp=WebGestaltR(enrichMethod="ORA",
                       organism=organism,
                       enrichDatabase=enrichDatabase,
                       interestGene=genes,
                       interestGeneType=interestGeneType,
                       referenceGene=referenceGene,
                       referenceGeneType=interestGeneType,
                       referenceSet=referenceSet,
                       isOutput=TRUE,
                       minNum=5,
                       outputDirectory=outputDirectory,
                       projectName=projectName,
                       sigMethod="top",
                       topThr=10)
    }
  }
}
webGestaltR_version<-paste0('WebGestaltR,v', packageVersion('WebGestaltR'))
writeLines(webGestaltR_version, 'WebGestaltR.version')

write.csv(allres, file=paste0(outFile, ".WebGestaltR.files.csv"), row.names=F)
