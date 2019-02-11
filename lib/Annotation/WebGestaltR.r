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

geneList<-read.table(geneFile,header=FALSE,sep="\t",stringsAsFactors=FALSE)
genes<-geneList$V1

enrichDatabases<-c("geneontology_Biological_Process", 
                   "geneontology_Cellular_Component", 
                   "geneontology_Molecular_Function",
                   "pathway_KEGG", 
                   "pathway_Wikipathway", 
                   "network_Kinase_target", 
                   "network_miRNA_target", 
                   "network_PPI_BIOGRID", 
                   "network_Transcription_Factor_target"
)

for(enrichDatabase in enrichDatabases){
  WebGestaltR(enrichMethod="ORA",organism=organism,
            enrichDatabase=enrichDatabase,interestGene=genes,
            interestGeneType=interestGeneType,referenceSet=referenceSet,
            isOutput=TRUE,
            outputDirectory=outputDirectory,projectName=paste0(sampleName, "_", enrichDatabase))
}
