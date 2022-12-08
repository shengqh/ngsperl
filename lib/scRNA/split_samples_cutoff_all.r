rm(list=ls()) 
outFile='P9112'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parFile1='/data/h_gelbard_lab/projects/20221129_9112_scRNA_human/hto_samples.txt'
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20221129_9112_scRNA_human/hto_samples_cutoff_all/result')

### Parameter setting end ###

source("split_samples_utils.r")
library(Seurat)
library(ggplot2)

#devtools::install_github("shengqh/cutoff")
#install.packages("bbmle")
library(zoo)
library(reshape2)
library(gridExtra)
library(ggExtra)

library(bbmle)
library(choisycutoff)

#still keep it here in case we need it later.
my_startval <- function(values,D1="normal",D2="normal",cutoff_point=0) {
  if(cutoff_point >0){
    thresh = cutoff_point
  }else{
    den <- tryCatch(
      expr = {
        density(values, bw="SJ")
      },
      error = function(e){ 
        density(values)
      }
    )
    w=1
    x=den$x
    y=den$y
    y.smooth=den$y
    n <- length(y.smooth)
    y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
    delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
    i.max <- which(delta <= 0) + w
    res=data.frame(x=x[i.max], i=i.max, y=y[i.max])
    res=res[res$x>0,]
    res=res[order(res$y, decreasing = T),]

    if(nrow(res)>2){
      res=res[1:2,]
    }

    xx=x[x>min(res$x) & x<max(res$x)]
    yy=y[x>min(res$x) & x<max(res$x)]
    yy.min=min(yy)
    ii.min=which(yy==yy.min)
    thresh=xx[ii.min]
  }

  sel <- values<thresh
  data1 <- values[sel]
  data2 <- values[!sel]
  lambda <- length(data1)/length(values)
  param1 <- MASS::fitdistr(data1,D1)$est
  param2 <- MASS::fitdistr(data2,D2)$est
  out <- c(param1,param2,lambda)
  names(out) <- c("mu1","sigma1","mu2","sigma2","lambda")
  return(out)
}

get_cutoff<-function(values, prefix, cutoff_startval=NA){
  my_out <- em(values,"normal","normal", cutoff_startval=cutoff_startval)
  cut_off <- cutoff(my_out)
  draw_cutoff(paste0(output_prefix, "_", tagname), values, cut_off[1])
  return(list(em_out=my_out, cut_off=cut_off))
}

files_lines=read.table(parSampleFile1, sep="\t")
files=split(files_lines$V1, files_lines$V2)

params_lines=read.table(parSampleFile3, sep="\t")
params=split(params_lines$V1, params_lines$V2)
params$hto_ignore_exists=ifelse(params$hto_ignore_exists=="0", FALSE, TRUE)
umap_min_dist=as.numeric(params$umap_min_dist)
umap_num_neighbors=as.numeric(params$umap_num_neighbors)

cutoffs=NULL
if(!is.na(params$cutoff_file)){
  if(file.exists(params$cutoff_file)){
    cutoffs=read.table(params$cutoff_file, sep="\t")
  }
}

has_hto_samples = FALSE
if(exists('parSampleFile4')){
  has_hto_samples=TRUE
  hto_samples_tbl = read.table(parSampleFile4, sep="\t", header=F)
  hto_samples = split(hto_samples_tbl$V1, hto_samples_tbl$V2)
}
if(!has_hto_samples & (parFile1 != "")){
  has_hto_samples=TRUE
  hto_samples_tbl = read.table(parFile1, sep="\t", header=T)
  hto_samples = split(hto_samples_tbl$Tagname, hto_samples_tbl$File)
}

idx=1
for(idx in c(1:length(files))){
  fname=names(files)[idx]
  output_prefix = paste0(fname, ".HTO")
  output_file=paste0(output_prefix, ".csv")
  
  #if(file.exists(output_file) & params$hto_ignore_exists){
  #  next
  #}

  rdsfile=files[[idx]]
  cat(fname, ":", rdsfile, " ...\n")

  if (has_hto_samples){
    cur_tags = hto_samples[[fname]]
  }else{
    cur_tags = NULL
  }

  obj=read_hto(rdsfile, output_prefix, cur_tags)
  
  tagnames=rownames(obj[["HTO"]])

  data <- FetchData(object=obj, vars=tagnames)
  write.csv(data, file=paste0(output_prefix, ".data.csv"))
  
  res=NULL
  tagname=tagnames[1]  
  for (tagname in tagnames) {
    values=data[,tagname]
    values=values[values>0]

    my_cutoff=lookup_cutoff(cutoffs, fname, tagname)
    if(my_cutoff == 0){
      cat(paste0("get cutoff of ", tagname, " ...\n"))
      my_cutoff=get_cutoff(values, paste0(output_prefix, "_", tagname))
      data[,paste0(tagname,"_pos")] = ifelse(data[,tagname]>my_cutoff$cut_off[1], tagname, "Negative")
      res<-rbind(res, data.frame(
        tagname=tagname,
        start.val=my_cutoff$em_out$start.val,
        mu1=my_cutoff$em_out$param["mu1"], 
        sigma1=my_cutoff$em_out$param["sigma1"], 
        mu2=my_cutoff$em_out$param["mu2"], 
        sigma2=my_cutoff$em_out$param["sigma2"], 
        cutoff=my_cutoff$cut_off[1]))
    }else{
      cat("use predefined cutoff=", my_cutoff, "\n")
      draw_cutoff(paste0(output_prefix, "_", tagname), values, my_cutoff)
      data[,paste0(tagname,"_pos")] = ifelse(data[,tagname]>my_cutoff, tagname, "Negative")
    }
  }
  if(!is.null(res)){
    write.csv(res, paste0(output_prefix, ".estimated.csv"))
  }
  
  cat(paste0("get classification ...\n"))
  class_names=paste0(tagnames, "_pos")
  data$HTO_classification=unlist(apply(data, 1, function(x){
    xx=unique(x[class_names])
    if (length(xx) > 1){
      xx = xx[xx != "Negative"]
      if (length(xx) > 1) {
        return("Doublet")
      }
    }
    return(xx)
  }))
  
  data$HTO_classification.global=unlist(apply(data, 1, function(x){
    xx=unique(x[class_names])
    if (length(xx) > 1){
      xx = xx[xx != "Negative"]
      if (length(xx) > 1) {
        return("Doublet")
      }else{
        return("Singlet")
      }
    }
    return("Negative")
  }))
  
  obj[["HTO_classification"]] = data$HTO_classification
  obj[["HTO_classification.global"]] = data$HTO_classification.global

  saveRDS(obj, file=paste0(output_prefix, ".umap.rds"))
  
  output_post_classification(obj, output_prefix, umap_min_dist=umap_min_dist, umap_num_neighbors=umap_num_neighbors, tagnames=tagnames)
}
