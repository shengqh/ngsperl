rm(list=ls()) 
outFile='P13988_rnaseq'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/data/cqs/references/gencode/GRCh38.p14/gencode.v48.annotation.gtf.map'
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/jennifer_pietenpol_projects/20251205_13988_Balko_hg38/20251205_13988_rnaseq/dexseq_count_table/result')

### Parameter setting end ###

library(data.table)
library(dplyr)

file_map=fread(parSampleFile1, data.table=FALSE, header=FALSE)
gene_tbl=fread(parFile1, data.table=FALSE, header=TRUE)
gene_map=split(gene_tbl$gene_name, gene_tbl$gene_id)

cat("Reading", nrow(file_map), "samples ...\n")

i=1
count_tbl=NULL
for(i in 1:nrow(file_map)){
  sample_name=file_map[i,2]
  count_file=file_map[i,1]
  cat("  Processing", sample_name, "from", count_file, "...\n")
  counts<-fread(file=count_file, sep="\t", data.table=FALSE, check.names=FALSE, quote="")
  colnames(counts)<-c("Feature", sample_name)

  # remove those rows: "_ambiguous" "_ambiguous_readpair_position" "_empty" "_lowaqual" "_notaligned"         
  counts = counts |> dplyr::filter(!grepl("^_", Feature))

  if(is.null(count_tbl)){
    count_tbl=counts
  }else{
    if(!all(count_tbl$Feature==counts$Feature)){
      stop("Exon IDs do not match for sample ", sample_name, "\n")
    }
    count_tbl=cbind(count_tbl, counts[,2,drop=FALSE])
  }
}

cat("Assigning gene information ...\n")

# Now assign gene information
count_tbl$Feature=gsub('"', '', count_tbl$Feature, fixed=TRUE)

exon_ids=do.call(rbind, strsplit(count_tbl$Feature, ":", fixed=TRUE)) |>
  as.data.frame(stringsAsFactors=FALSE) |>
  dplyr::rename(Feature_gene_ids=V1, Feature_exon_num=V2)

count_tbl=cbind(exon_ids, count_tbl)

multi_gene_ids = count_tbl |> dplyr::filter(grepl("\\+", Feature_gene_ids)) |>
  dplyr::select(Feature_gene_ids) |>
  unique() |>
  dplyr::pull()

multi_genes=lapply(multi_gene_ids, function(x){
  parts=unlist(strsplit(x, "+", fixed=TRUE))
  gene_names=gene_map[parts]
  unique_gene_names=unique(unlist(gene_names))
  if(length(unique_gene_names)==1){
    return(unique_gene_names)
  }else{
    return(paste(gene_names, collapse="+"))
  }
}) |> unlist()

multi_gene_map=setNames(multi_genes, multi_gene_ids)
all_gene_map=c(gene_map, multi_gene_map)

stopifnot(all(count_tbl$Feature_gene_ids %in% names(all_gene_map)))

count_tbl$Feature_gene_names=unlist(all_gene_map[count_tbl$Feature_gene_ids])

count_tbl = count_tbl |>
  dplyr::select(Feature, Feature_gene_ids, Feature_gene_names, Feature_exon_num, dplyr::everything())

count_file=paste0(outFile, ".exon.count")
cat("Writing count table to", count_file, "...\n") 
write.table(count_tbl, file=count_file, quote=FALSE, sep="\t", row.names=FALSE)
