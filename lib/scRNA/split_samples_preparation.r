
library(Seurat)

params_lines=read.table(parSampleFile2, sep="\t")
params=split(params_lines$V1, params_lines$V2)
params$hto_non_zero_percentage=as.numeric(params$hto_non_zero_percentage)
params$nFeature_cutoff_min=as.numeric(params$nFeature_cutoff_min)

files_lines=read.table(parSampleFile1, sep="\t")
files=split(files_lines$V1, files_lines$V2)

cutoff_tbl<-NULL

output_tag_dist<-function(htos, filename){
  obj <- CreateSeuratObject(counts = htos, assay="HTO")
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
  DefaultAssay(object = obj) <- "HTO"
  tagnames=rownames(obj[["HTO"]])
  n_col=ceiling(sqrt(length(tagnames)))
  n_row=ceiling(length(tagnames) / n_col)
  width=max(3200, n_col * 3000 + 200)
  height=max(3000, n_row * 3000)
  png(filename, width=width, height=height, res=300)
  rplot(object=obj, assay="HTO", features = tagnames, identName="orig.ident", n_row=n_row)
  dev.off()
}

cname=names(files)[1]
for(cname in names(files)){
  cfiles=files[[cname]]
  if(length(cfiles) == 1){
    cat("preparing", cname, ":", cfiles, " ...\n")
    if(grepl('.h5$', cfiles)){
      sdata<-Read10X_h5(cfiles)
      exp<-sdata[[1]]
      writeLines(colnames(exp),con = paste0(cname, ".barcodes.tsv"))
      htos<-sdata[[2]]
      htos<-as.matrix(htos)
    }else{
      stop(paste0("Unknown file format:", cfiles))
    }
  }else{
    hfile=cfiles[[1]]
    cat("preparing", cname, ":", hfile, " ...\n")
    efile=cfiles[[2]]
    if(grepl('hto_mtx.rds$', hfile)){
      htos<-readRDS(hfile)
      htos<-as.matrix(htos)
      #tag number should less than cell number
      if(ncol(htos) < nrow(htos)){
        htos=t(htos)
      }
      exp<-readRDS(efile)
      cells.use <- intersect(colnames(exp), colnames(htos))
      exp<-exp[, cells.use]
      htos<-htos[, cells.use]
    }else{
      stop(paste0("Unknown file format:", cfiles))
    }
  }

  write.csv(htos, paste0(cname, ".raw.alltags.csv"), row.names=T)

  if (!is.na(params$hto_regex) & params$hto_regex != "" ) {
    htos<-htos[grepl(params$hto_regex, rownames(htos)),]
    if (nrow(htos) == 0){
      stop(paste0("Cannot find hashtag based on regex ", params$hto_regex, " for tags ", paste(rownames(mat), collapse=",")))
    }
    cat("After hash tag regex filter: ", paste(rownames(htos), collapse=","), "\n")
  }

  output_tag_dist(htos, paste0(cname, ".raw.alltags.png"))

  hto.exp <- CreateSeuratObject(counts = exp, min.features = params$nFeature_cutoff_min)
  cells.valid<-colnames(hto.exp)
  htos<-htos[,cells.valid]
  
  write.csv(htos, paste0(cname, ".alltags.csv"), row.names=T)

  if(params$hto_non_zero_percentage > 0){
    rowsum<-apply(htos>0, 1, sum)
    htos<-htos[rowsum > (ncol(htos) * params$hto_non_zero_percentage),]
    cat("After zero count filter: ", paste(rownames(htos), collapse=","), "\n")
  }
  
  rownames(htos)<-gsub("^TotalSeqC_", "", rownames(htos))
  rownames(htos)<-gsub("^TotalSeq_", "", rownames(htos))
  rownames(htos)<-gsub('.TotalSeqC$', "", rownames(htos))
  
  cat("After name clean: ", paste(rownames(htos), collapse=","), "\n")
  
  obj <- CreateSeuratObject(counts = htos, assay="HTO")
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
  DefaultAssay(object = obj) <- "HTO"
  
  tagnames=rownames(obj[["HTO"]])
  
  cutoff_tbl<-rbind(cutoff_tbl, data.frame("cutoff"=0, "tagname"=tagnames, "filename"=cname))
  
  n_col=ceiling(sqrt(length(tagnames)))
  n_row=ceiling(length(tagnames) / n_col)

  width=max(3200, n_col * 3000 + 200)
  height=max(3000, n_row * 3000)
  png(paste0(cname, ".tag.dist.png"), width=width, height=height, res=300)
  rplot(object=obj, assay="HTO", features = tagnames, identName="orig.ident", n_row=n_row)
  dev.off()

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

  saveRDS(htos, paste0(cname, ".hto.rds"))
}

write.table(cutoff_tbl, paste0(outFile, ".cutoff_template.txt"), row.names = F, col.names = F, sep="\t", quote=F)
