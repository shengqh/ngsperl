options(bitmapType='cairo')

library(WebGestaltR)

args = commandArgs(trailingOnly=TRUE)
organism = args[1] #hsapiens
sampleName=args[2]
geneFile = args[3]
outputDirectory = args[4]
interestGeneType = args[5]
referenceSet = args[6]

cat("organism=", organism, "\n")
cat("sampleName=", sampleName, "\n")
cat("geneFile=", geneFile, "\n")
cat("outputDirectory=", outputDirectory, "\n")

if(!exists("interestGeneType")){
  interestGeneType="genesymbol"
}

if(!exists("referenceSet")){
  referenceSet="genome"
}
cat("interestGeneType=", interestGeneType, "\n")
cat("referenceSet=", referenceSet, "\n")

info = file.info(geneFile)
if (all(is.na(info$size))) {
  stop(paste0("Gene file is not exist: ", geneFile, "\n"))
}

if (all(info$size == 0)) {
  stop(paste0("Gene file is empty: ", geneFile, "\n"))
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


enrichDatabases<-c("geneontology_Biological_Process", 
                   "geneontology_Cellular_Component", 
                   "geneontology_Molecular_Function",
                   "pathway_KEGG" 
                   #"pathway_Wikipathway", 
                   #"network_miRNA_target",
                   #"network_PPI_BIOGRID", 
                   #"network_Transcription_Factor_target"
)

for(enrichDatabase in enrichDatabases){
  temp=WebGestaltR(enrichMethod="ORA",organism=organism,
            enrichDatabase=enrichDatabase,interestGene=genes,
            interestGeneType=interestGeneType,referenceSet=referenceSet,
            isOutput=TRUE,minNum=5,
            outputDirectory=outputDirectory,projectName=paste0(sampleName, "_", enrichDatabase))
  if (is.null(temp)) { #no enrichment. report top 5 categories
    warning(paste0("No significant category (FDR<=0.05) identified in ",enrichDatabase,", reporting top 10 categories instead."))
    temp=WebGestaltR(enrichMethod="ORA",organism=organism,
                     enrichDatabase=enrichDatabase,interestGene=genes,
                     interestGeneType=interestGeneType,referenceSet=referenceSet,
                     isOutput=TRUE,minNum=5,
                     outputDirectory=outputDirectory,projectName=paste0(sampleName, "_", enrichDatabase),
                     sigMethod="top",topThr=10)
  }
}

webGestaltR_version<-paste0('WebGestaltR,v', packageVersion('WebGestaltR'))
writeLines(webGestaltR_version, 'WebGestaltR.version')
