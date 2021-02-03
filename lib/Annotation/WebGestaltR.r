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
  geneList<-read.csv(geneFile,header=FALSE,stringsAsFactors=FALSE)
} else {
  geneList<-read.table(geneFile,header=FALSE,sep="\t",stringsAsFactors=FALSE)
}

if (ncol(geneList)==1) {
  genes<-geneList$V1
} else { #try to find gene column
  geneInd=grep("Gene|gene",geneList[1,])
  if (length(geneInd)>0) {
    genes<-geneList[,geneInd[1]]
    genes<-genes[2:length(genes)]
  } else { #guess gene column, by contents with both number and character (so not all numeric data, can be gene IDs)
    geneInd=which.max(apply(geneList[-1,],2,function(x) length(intersect(grep("[a-zA-Z][a-zA-Z]",x),grep("[0-9]",x)))))
    genes<-geneList[,geneInd]
  }
}

# if(grepl("Gene", genes[1])){
#   genes<-genes[2:length(genes)]
# }

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
  WebGestaltR(enrichMethod="ORA",organism=organism,
            enrichDatabase=enrichDatabase,interestGene=genes,
            interestGeneType=interestGeneType,referenceSet=referenceSet,
            isOutput=TRUE,
            outputDirectory=outputDirectory,projectName=paste0(sampleName, "_", enrichDatabase))
}

webGestaltR_version<-paste0('WebGestaltR,v', packageVersion('WebGestaltR'))
writeLines(webGestaltR_version, 'WebGestaltR.version')
