options(bitmapType='cairo')

library(WebGestaltR)

args = commandArgs(trailingOnly=TRUE)
sampleName=args[1]
geneFile = args[2]
outputDirectory = args[3]

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
  WebGestaltR(enrichMethod="ORA",organism="hsapiens",
            enrichDatabase=enrichDatabase,interestGene=genes,
            interestGeneType="genesymbol",referenceSet="genome",
            is.output=TRUE,
            outputDirectory=outputDirectory,projectName=paste0(sampleName, "_", enrichDatabase))
}
