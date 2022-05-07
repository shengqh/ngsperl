
library(Seurat)
library(ggplot2)

#devtools::install_github("shengqh/cutoff")
#install.packages("bbmle")
library(zoo)
library(reshape2)
library(gridExtra)
library(ggExtra)

library(autothresholdr)

#https://stats.stackexchange.com/questions/500148/difficulty-splitting-bimodal-data
find_cutoff<-function(values){
  ivalues=floor(values * 100)
  cutoff<-autothresholdr::auto_thresh(ivalues, "InterModes")
  return(cutoff/100)
}

draw_cutoff<-function(prefix, values, cut_off){
  png(paste0(prefix, ".cutoff.png"), width=2000, height=1600, res=300)
  his<-hist(values,200,F,xlab="concentration",ylab="density", main=NULL,col="grey")
  lines(density(values),lwd=1.5,col="blue")
  abline(v=cut_off,lwd=1.5,col="red")

  minb=min(his$breaks)
  maxb=max(his$breaks)
  x=minb + (maxb-minb) * 3 /4
  
  y=max(his$density)
  legend(x=x, y=y, legend=c("density", "cutoff"), col=c("blue", "red"), lty=c(1,2))
  dev.off()
}

files_lines=read.table(parSampleFile1, sep="\t")
files=split(files_lines$V1, files_lines$V2)

params_lines=read.table(parSampleFile3, sep="\t")
params=split(params_lines$V1, params_lines$V2)
params$hto_ignore_exists=ifelse(params$hto_ignore_exists=="0", FALSE, TRUE)
umap_min_dist=as.numeric(params$umap_min_dist)

if(file.exists(params$cutoff_file)){
  cutoffs=read.table(params$cutoff_file, sep="\t")
}else{
  cutoffs=NULL
}

idx=14
for(idx in c(1:length(files))){
  fname=names(files)[idx]
  output_prefix = paste0(fname, ".HTO")
  output_file=paste0(output_prefix, ".csv")
  
  #if(file.exists(output_file) & params$hto_ignore_exists){
  #  next
  #}

  rdsfile=files[[idx]]
  cat(fname, ":", rdsfile, " ...\n")

  obj=read_hto(rdsfile, output_prefix)
  
  tagnames=rownames(obj[["HTO"]])
  data <- FetchData(object=obj, vars=tagnames)
  write.csv(data, file=paste0(output_prefix, ".data.csv"))
  
  res=NULL
  tagname=tagnames[1]  
  for (tagname in tagnames) {
    values=data[,tagname]
    values=values[values>0]

    my_cutoff=lookup_cutoff(cutoffs, fname, tagname)

    cat(paste0("get cutoff of ", tagname, " ...\n"))
    if(my_cutoff == 0){
      my_cutoff=find_cutoff(values)
      cat("find cutoff=", my_cutoff, "\n")
    }else{
      cat("use predefined cutoff=", my_cutoff, "\n")
    }
    draw_cutoff(paste0(output_prefix, "_", tagname), values, my_cutoff)
    
    data[,paste0(tagname,"_pos")] = ifelse(data[,tagname]>my_cutoff, tagname, "Negative")
    res<-rbind(res, data.frame(
      sample=fname,
      tag=tagname,
      cutoff=my_cutoff))
  }
  write.csv(res, paste0(output_prefix, ".estimated.csv"), row.names=F)
  
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
  
  output_post_classification(obj, output_prefix, umap_min_dist=umap_min_dist)
}
