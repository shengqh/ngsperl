rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/shengq2/20230115_combined_scRNA_hg38/hto_samples_preparation/result')

### Parameter setting end ###

source("split_samples_utils.r")
library('Seurat')
library('R.utils')
library('reshape2')
library('Matrix')

params_lines=read.table(parSampleFile2, sep="\t")
params=split(params_lines$V1, params_lines$V2)
params$hto_non_zero_percentage=as.numeric(params$hto_non_zero_percentage)
params$nFeature_cutoff_min=as.numeric(params$nFeature_cutoff_min)
params$hto_filter_by_exp=is_one(params$hto_filter_by_exp)
params$hto_min_count=as.numeric(params$hto_min_count)

files_lines=read.table(parSampleFile1, sep="\t")
files=split(files_lines$V1, files_lines$V2)

has_raw_files = parSampleFile3 != ""
if(has_raw_files){
  files_lines=read.table(parSampleFile3, sep="\t")
  raw_files=split(files_lines$V1, files_lines$V2)
}
cutoff_tbl<-NULL

cname = names(files)[1]
for(cname in names(files)){
  cfiles = files[[cname]]
  res_lst = read_hto_file(cname, cfiles)
  htos=res_lst$htos
  exp=res_lst$exp
  rm(res_lst)
  write.csv(htos, paste0(cname, ".alltags_raw.csv"), row.names=T)

  if(has_raw_files){
    raw_lst = read_hto_file(cname, raw_files[[cname]])
    htos=raw_lst$htos
    rm(raw_lst)
  }

  if (!is.na(params$hto_regex) & params$hto_regex != "" ) {
    htos<-htos[grepl(params$hto_regex, rownames(htos)),]
    if (nrow(htos) == 0){
      stop(paste0("Cannot find hashtag based on regex ", params$hto_regex, " for tags ", paste(rownames(mat), collapse=",")))
    }
    cat("After hash tag regex filter: ", paste(rownames(htos), collapse=","), "\n")
  }

  obj <- CreateSeuratObject(counts = htos, assay="HTO")

  htos <- obj@assays$HTO@counts

  cat("Final tagnames:", paste0(rownames(htos), collapse = ", "), "\n")
  
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
  DefaultAssay(object = obj) <- "HTO"

  output_tag_dist(obj, paste0(cname, ".alltags_raw.png"))

  if(params$hto_filter_by_exp){
    hto.exp <- CreateSeuratObject(counts = exp, min.features = params$nFeature_cutoff_min)
    cells.valid<-colnames(hto.exp)
    htos<-htos[,cells.valid]
  }else{
    htos<-htos[,colSums(htos) >= params$hto_min_count,drop=F]
  }
  
  write.csv(htos, paste0(cname, ".alltags.csv"), row.names=T)

  if(params$hto_non_zero_percentage > 0){
    rowsum<-apply(htos>0, 1, sum)
    htos<-htos[rowsum > (ncol(htos) * params$hto_non_zero_percentage),]

    if(nrow(htos) == 0){
      cat("After zero count filter, no tag left, ignored\n")
      next
    }

    cat("After zero count filter: ", paste(rownames(htos), collapse=","), "\n")
    obj<-subset(obj, features=rownames(htos))
  }
  
  obj$filtered<-colnames(obj) %in% colnames(htos)
  
  tagnames=rownames(obj[["HTO"]])
  
  cutoff_tbl<-rbind(cutoff_tbl, data.frame("cutoff"=0, "tagname"=tagnames, "filename"=cname))

  n_col=ceiling(sqrt(length(tagnames)))
  n_row=ceiling(length(tagnames) / n_col)
  if(has_raw_files){
    identName='filtered'
  }else{
    identName='orig.ident'
  }

  width=min(20000, n_col * 3000 + 200)
  height=min(20000, n_row * 1500 * ifelse(has_raw_files, 2, 1))
  png(paste0(cname, ".tag.dist.png"), width=width, height=height, res=300)
  rplot(object=obj, assay="HTO", features = tagnames, identName=identName, n_row=n_row)
  dev.off()

  obj<-subset(obj, cells=colnames(htos))

  if (length(tagnames) == 2) {
    png(paste0(cname, ".tag.point.png"), width=2000, height=1800, res=300)
    print(FeatureScatter(object = obj, feature1 = tagnames[1], feature2 = tagnames[2], cols = "black"))
    dev.off()
  }
  # fdata<-FetchData(obj, assay="HTO", features=tagnames)
  # mt<-data.frame(UMAP_1=obj@reductions$umap@cell.embeddings[,1], 
  #              UMAP_2=obj@reductions$umap@cell.embeddings[,2],
  #              Sample=obj$orig.ident,
  #              batch=obj$batch)

  saveRDS(obj, paste0(cname, ".hto.rds"))

  writeLines(colnames(obj), paste0(cname, ".barcodes.tsv"))

  counts=as.sparse(htos)

  save_to_matrix(counts, cname);
}

write.table(cutoff_tbl, paste0(outFile, ".cutoff_template.txt"), row.names = F, col.names = F, sep="\t", quote=F)
