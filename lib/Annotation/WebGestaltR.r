options(bitmapType='cairo')

library(WebGestaltR)

args = commandArgs(trailingOnly=TRUE)
organism = args[1] #hsapiens
sampleName=args[2]
geneFile = args[3]
outputDirectory = args[4]

cat("organism=", organism, "\n")
cat("sampleName=", sampleName, "\n")
cat("geneFile=", geneFile, "\n")
cat("outputDirectory=", outputDirectory, "\n")

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
            interestGeneType="genesymbol",referenceSet="genome",
            is.output=TRUE,
            outputDirectory=outputDirectory,projectName=paste0(sampleName, "_", enrichDatabase))
}
