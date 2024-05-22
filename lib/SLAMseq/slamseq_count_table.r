rm(list=ls()) 
outFile='P11389_SLAMseq'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/data/cqs/references/gencode/GRCm38.p6/gencode.vM25.annotation.gtf'
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20240520_11389_SLAMseq_mm10/T11_count_table/result')

### Parameter setting end ###

source("countTableVisFunctions.R")
library(data.table)
library(DESeq2)

files=read_file_map(parSampleFile1)

gene_map = fread(parFile1, data.table=FALSE)
id_tbl = gene_map |> 
  dplyr::filter(V3=="transcript") |> 
  dplyr::select(V9) |>
  dplyr::mutate(
    g1 = gsub('"; transcript_id.+', '', V9),
    gene_id = gsub('gene_id "', '', g1),
    g2 = gsub('"; gene_type.+', '', V9),
    transcript_id = gsub('.+"', '', g2),
    g3 = gsub('"; transcript_type.+', '', V9),
    gene_name = gsub('.+"', '', g3)
  ) |>
  dplyr::select(gene_id, transcript_id, gene_name)

all_read_df=NULL
tc_read_df=NULL
sample = names(files)[1]
for(sample in names(files)){
  df = fread(files[[sample]])
  df$gene_name = gsub("_3utr","",df$gene_name)

  read_count_df = df[,c("gene_name","readCount")]
  colnames(read_count_df) = c("transcript_id",sample)
  if(is.null(all_read_df)){
    all_read_df = read_count_df
  }else{
    all_read_df = merge(all_read_df, read_count_df, by="transcript_id", all=TRUE)
  }

  tcread_count_df = df[,c("gene_name","tcReadCount")]
  colnames(tcread_count_df) = c("transcript_id",sample)
  if(is.null(tc_read_df)){
    tc_read_df = tcread_count_df
  }else{
    tc_read_df = merge(tc_read_df, tcread_count_df, by="transcript_id", all=TRUE)
  }  
}

stopifnot(all(all_read_df$transcript_id %in% id_tbl$transcript_id))

all_read_df=all_read_df[rowSums(all_read_df |> dplyr::select(-transcript_id)) > 0,]
all_read_tbl = merge(id_tbl, all_read_df, by="transcript_id", all.x=FALSE, all.y=TRUE) |>
  dplyr::select(-transcript_id)
all_read_tbl = aggregate(. ~ gene_id + gene_name, all_read_tbl, sum, na.rm=TRUE) |>
  dplyr::rename(Feature=gene_id, Feature_gene_name=gene_name)

write.csv(all_read_tbl, file=paste0(outFile, ".all_read.csv"), quote=FALSE, row.names=FALSE)

tc_read_df=tc_read_df[rowSums(tc_read_df |> dplyr::select(-transcript_id)) > 0,]
tc_read_tbl = merge(id_tbl, tc_read_df, by="transcript_id", all.x=FALSE, all.y=TRUE) |>
  dplyr::select(-transcript_id)
tc_read_tbl = aggregate(. ~ gene_id + gene_name, tc_read_tbl, sum, na.rm=TRUE) |>
  dplyr::rename(Feature=gene_id, Feature_gene_name=gene_name)
write.csv(tc_read_tbl, file=paste0(outFile, ".tc_read.csv"), quote=FALSE, row.names=FALSE)

myEstimateSizeFactors<-function(dds){
  if(exists("librarySize")){
    cat("Estimate size factor based on library size\n")
    curLibrarySize<-librarySize[colnames(dds)]
    #based on DESeq2 introduction
    curSizeFactor<- curLibrarySize / exp(mean(log(curLibrarySize)))
    sizeFactors(dds)<-curSizeFactor
  }else{
    cat("Estimate size factor based on reads\n")
    sfres<-try(dds<-estimateSizeFactors(dds))
    if (class(sfres) == "try-error") {
      library(edgeR)
      countNum<-counts(dds)
      y<-calcNormFactors(countNum, methold="TMM")
      cs<-colSums(countNum)
      cs<-cs / median(cs)
      sf<-y * cs
      sizeFactors(dds)<-sf
    }
  }
  return(dds)
}

countNum=all_read_df |> tibble::column_to_rownames("transcript_id") |> as.matrix()
designData=data.frame(Sample=colnames(countNum))
dds=DESeqDataSetFromMatrix(countData = countNum,
                          colData = designData,
                          design = ~1) 
dds=myEstimateSizeFactors(dds)
sf=data.frame(Sample=colnames(countNum), sizeFactor=sizeFactors(dds))
write.csv(sf, file=paste0(outFile, ".sizeFactor.csv"), quote=FALSE, row.names=FALSE)

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
