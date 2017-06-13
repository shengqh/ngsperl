options(bitmapType='cairo')

#inputfile<-"e:/temp/genename.txt"
#outputfile<-"e:/temp/genelocation.bed"
host<-"useast.ensembl.org"
dataset<-"hsapiens_gene_ensembl"

resultFile<-outFile
inputFile<-parFile1

library(biomaRt)

genetable<-read.table(inputfile, sep="\t", header=F, check.names = F, row.names=1)
genesymbols<-rownames(genetable)

gss<-sub("_.*$", "", genesymbols, perl=TRUE) 
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host=host, dataset=dataset)
genepos <- lapply(gss, function(x){
  getBM(attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
        filters="external_gene_name",
        values = x,
        mart=mart)})
genepos <- do.call("rbind", genepos)
colnames(genepos)<-c("chr",	"s1",	"s2", "geneid")
write.table(genepos, file=outputfile, row.names = F,col.names = F,  quote = F, sep="\t")
