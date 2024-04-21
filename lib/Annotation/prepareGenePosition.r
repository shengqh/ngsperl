options(bitmapType='cairo')

library(biomaRt)
library(gtools)

args = commandArgs(trailingOnly=TRUE)
datasetName = args[1]
outputFile = args[2]
isGFF = args[3] == "1"
addChr = args[4] == "1"
genes=args[5:length(args)]

#datasetName="hg19"
#isGFF=1
#outputFile="h:/temp.pos"
#genes=c("KDR")

datasets=c("hsapiens_gene_ensembl")
names(datasets)=c(datasetName)

symbolColumns=c("hgnc_symbol")
names(symbolColumns)=names(datasets)

if(datasetName %in% names(datasets)){
  host<-"www.ensembl.org"
  dataset<-datasets[datasetName]
  symbolColumn<-symbolColumns[datasetName]
  
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host=host, dataset=dataset)
  allgenepos <- lapply(genes, function(x){
    getBM(attributes = c("chromosome_name", "start_position", "end_position", symbolColumn, "strand"),
          filters=symbolColumn,
          values = x,
          mart=mart,
          useCache=FALSE)})
  allgenepos <- do.call("rbind", allgenepos)
  rownames(allgenepos)<-allgenepos$ID
  colnames(allgenepos)<-c("chr",  "start",  "end", "gene", "strand")
  allgenepos$strand[allgenepos$strand==1] = '+'
  allgenepos$strand[allgenepos$strand==-1] = '-'
  
  if(addChr){
    allgenepos$chr = paste0("chr", allgenepos$chr)
  }
  
  ll = unlist(lapply(allgenepos$chr, nchar))
  allgenepos = allgenepos[ll < 6,]

  if(isGFF){
    allgenepos$unknown1="GENE"
    allgenepos$unknown2='.'
    allgenepos$unknown3='.'
    allgenepos=allgenepos[,c(1,4,6,2,3,7,5,8)]
    write.table(allgenepos, file=outputFile, row.names = F,col.names = F,  quote = F, sep="\t")
  }else{
    write.table(allgenepos, file=outputFile, row.names = F,col.names = T,  quote = F, sep="\t")
  }
}else{
  stop(paste0("dataset ", datasetName, " is not in allowed list ", names(datasets)))
}
