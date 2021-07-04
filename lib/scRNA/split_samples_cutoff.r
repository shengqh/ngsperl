library(Seurat)
library(ggplot2)

#devtools::install_github("shengqh/cutoff")
#install.packages("bbmle")
library('choisycutoff')
library(zoo)
library(reshape2)
library(gridExtra)
library(ggExtra)

rplot<-function(object, features, assay, identName, withAllCells=FALSE){
  DefaultAssay(object = object) <- assay
  data <- FetchData(object = object, vars = c(features, identName))
  mdata<-melt(data, id.vars=identName)
  if (withAllCells) {
    mdata2<-mdata
    mdata2[,1] = "All cells"
    mdata<-rbind(mdata, mdata2)
  }
  gfinal=list()
  for(feature in features){
    ddata=mdata[mdata$variable==feature,]
    g<-ggplot(ddata, aes_string(x="value")) + 
      geom_histogram(aes(y=..density..), bins=50, colour="black", fill="white", position="identity") + 
      geom_density(color="red") +
      facet_grid(reformulate(".", identName), scale="free_y") + 
      xlab(feature) + theme_bw() + theme(strip.background=element_rect(colour="black", fill=NA),
                                         strip.text = element_text(size = 24),
                                         axis.text=element_text(size=18),
                                         axis.title=element_text(size=24))
    if (feature != features[1]){  
      g = g + ylab("")
    }
    gfinal = append(gfinal, list(g))
  }
  grid.arrange(grobs=gfinal, nrow=1)
}

my_startval <- function(values,D1="normal",D2="normal") {
  den=density(values)
  
  w=1
  x=den$x
  y=den$y
  y.smooth=den$y
  n <- length(y.smooth)
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  res=data.frame(x=x[i.max], i=i.max, y=y[i.max])
  res=res[res$x > 0,]
  res=res[order(res$y, decreasing = T),]

  #the highest peaks should be the negative one, positive one is always in the right side.
  res=res[res$x>=res$x[1],]

  if(nrow(res)>2){
    res=res[1:2,]
  }
  res=res[order(res$x),]
  
  xx=x[x>min(res$x) & x<max(res$x)]
  yy=y[x>min(res$x) & x<max(res$x)]
  yy.min=min(yy)
  ii.min=which(yy==yy.min)
  thresh=xx[ii.min]
  
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

t=1e-64
my_em<-function(values, data_name="em", D1="normal", D2="normal", t=1e-64){
  start <- as.list(my_startval(values, D1, D2))
  #print(start)
  
  D1b <- choisycutoff:::hash[[D1]]
  D2b <- choisycutoff:::hash[[D2]]
  lambda0 <- 0
  with(start, {
    while (abs(lambda0 - mean(lambda)) > t) {
      lambda <- mean(lambda)
      lambda0 <- lambda
      distr1 <- lambda * D1b(values, mu1, sigma1)
      distr2 <- (1 - lambda) * D2b(values, mu2, sigma2)
      lambda <- distr1/(distr1 + distr2)
      mLL2 <- function(mu1, sigma1, mu2, sigma2) return(choisycutoff:::mLL(mu1, 
                                                                     sigma1, mu2, sigma2, lambda, values, D1b, D2b))
      start <- as.list(log(c(mu1 = mu1, sigma1 = sigma1, 
                             mu2 = mu2, sigma2 = sigma2)))
      out <- bbmle::mle2(mLL2, start, "Nelder-Mead")
      coef <- out@coef
      coef_n <- names(coef)
      names(coef) <- NULL
      for (i in 1:4) assign(coef_n[i], exp(coef[i]))
    }
    out <- list(lambda = lambda, param = exp(out@coef), D1 = D1, 
                D2 = D2, deviance = out@min, data = values, data_name = data_name, 
                out = out, t = t)
    class(out) <- "em"
    return(out)
  })
}

my_cutoff<-function (object, t = 1e-64, nb = 10, distr = 2, type1 = 0.05, level = 0.95) 
{
  coef <- object$out@coef
  the_names <- names(coef)
  coef <- exp(mc2d::rmultinormal(nb, coef, as.vector(object$out@vcov)))
  coef <- as.list(data.frame(t(coef)))
  coef <- lapply(coef, function(x) {
      names(x) <- the_names
      return(as.list(x))
  })
  out <- sapply(coef, function(x) choisycutoff:::lci0(x, mean(object$lambda), 
      choisycutoff:::hash[[object$D1]], choisycutoff:::hash[[object$D2]], object$data, t))
  lambda <- rnorm(nb, out[1, ], out[2, ])
  coef <- sapply(coef, function(x) unlist(x))
  the_names <- c(rownames(coef), "lambda")
  coef <- rbind(coef, lambda)
  coef <- lapply(as.data.frame(coef), function(x) {
      names(x) <- the_names
      return(as.list(x))
  })
  out <- sapply(coef, function(x) with(x, choisycutoff:::cutoff0(mu1, 
      sigma1, mu2, sigma2, lambda, object$D1, object$D2, distr, type1)))
  out <- MASS::fitdistr(out, "normal")
  the_mean <- out$estimate["mean"]
  level <- (1 - level)/2
  level <- c(level, 1 - level)
  ci <- the_mean + qt(level, Inf) * out$sd["mean"]
  out <- c(the_mean, ci)
  names(out) <- c("Estimate", paste(100 * level, "%"))
  return(out)
}

get_cutoff<-function(values, prefix){
  my_out <- my_em(values,"normal","normal")
  cut_off <- my_cutoff(my_out)
  png(paste0(prefix, ".cutoff.png"), width=2000, height=1600, res=300)
  hist(values,200,F,xlab="concentration",ylab="density", main=NULL,col="grey")
  lines(density(values),lwd=1.5,col="blue")
  lines(my_out,lwd=1.5,col="red")
  abline(v=cut_off[1],lwd=1.5,col="brown")
  dev.off()
  return(cut_off[1])
}

split<-function(h5file, output_prefix, hashtag_regex=NA) {
  if(grepl('.rds$', h5file)){
    meta<-readRDS(h5file)
    meta<-as.matrix(meta)
    #tag number should less than cell number
    if(ncol(meta) < nrow(meta)){
      meta=t(as.matrix(meta))
    }
  }else{
    sdata<-Read10X_h5(h5file)
    meta<-sdata[[2]]
  }
  mat<-as.matrix(meta)
  rowsum<-apply(mat>0, 1, sum)
  mat<-mat[rowsum > (ncol(mat) / 2),]
  write.csv(mat, file=paste0(output_prefix, ".alltags.exp.csv"))
  
  if (!is.na(hashtag_regex)) {
    htos<-mat[grepl(hashtag_regex, rownames(mat)),]
    if (nrow(htos) == 0){
      stop(paste0("Cannot find hashtag based on regex ", hashtag_regex, " for tags ", paste(rownames(mat), collapse=",")))
    }
  }else{
    htos<-mat
  }
  rownames(htos)<-gsub("^TotalSeqC_", "", rownames(htos))
  rownames(htos)<-gsub("^TotalSeq_", "", rownames(htos))
  rownames(htos)<-gsub('.TotalSeqC$', "", rownames(htos))

  write.csv(htos, file=paste0(output_prefix, ".hto.exp.csv"))
  
  pbmc.hashtag <- CreateSeuratObject(counts = htos)
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  pbmc.hashtag <- NormalizeData(pbmc.hashtag, normalization.method = "CLR")
  
  tagnames=rownames(pbmc.hashtag)
  
  width=max(10, length(tagnames) * 5)
  pdf(paste0(output_prefix, ".tag.dist.pdf"), width=width, height=6)
  rplot(pbmc.hashtag, assay="RNA", features = tagnames, identName="orig.ident")
  dev.off()
  
  data <- FetchData(object=pbmc.hashtag, vars=tagnames)
  write.csv(data, file=paste0(output_prefix, ".hto.norm_exp.csv"))

  tagname=tagnames[1]  
  for (tagname in tagnames) {
    values=data[,tagname]
    cat(paste0("get cutoff of ", tagname, " ...\n"))
    cutoff=get_cutoff(values, paste0(output_prefix, "_", tagname))
    data[,paste0(tagname,"_pos")] = ifelse(data[,tagname]>cutoff, tagname, "Negative")
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
  
  pbmc.hashtag[["HTO_classification"]] = data$HTO_classification
  pbmc.hashtag[["HTO_classification.global"]] = data$HTO_classification.global
  
  width=max(1600, length(tagnames) * 600)
  height=1400

  png(paste0(output_prefix, ".class.dist.png"), width=width, height=height)
  rplot(pbmc.hashtag, assay = "RNA", features = tagnames, identName="HTO_classification")
  dev.off()
  
  if (length(tagnames) == 2) {
    pdf(paste0(output_prefix, ".class.point.png"), width=width, height=height)
    print(FeatureScatter(object = pbmc.hashtag, feature1 = tagnames[1], feature2 = tagnames[2],group.by="HTO_classification"))
    dev.off()
  }
  
  tmat=data.frame(t(mat))
  tmat$HTO = pbmc.hashtag$HTO_classification
  tmat$HTO.global = pbmc.hashtag$HTO_classification.global
  tmat$nCount_RNA = pbmc.hashtag$nCount_RNA
  tmat$nFeature_RNA = pbmc.hashtag$nFeature_RNA
  write.csv(tmat, paste0(output_prefix, ".csv"))
  
  width=length(tagnames) * 3
  pdf(paste0(output_prefix, ".class.RNA_Feature.pdf"), width=width, height=4)
  g<-ggplot(tmat, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point() + facet_grid(~HTO.global) + theme_bw() + theme(strip.background = element_blank())
  print(g)
  dev.off()
}

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  h5file = "C:/projects/data/cqs/alexander_gelbard_data/AG_5126_10X/Count/5126-AG-3/filtered_feature_bc_matrix.h5"
  output_prefix = "C:/projects/scratch/cqs/shengq2/papers/20210703_scrna_hto/hto_samples_cutoff/result/LTS4/LTS4.HTO"
  hashtag_regex="Hashtag|TotalSeqC"
  #h5file = "C:/Users/sheng/projects/paula_hurley/20201208_scRNA_split/filtered_feature_bc_matrix.h5"
  #output_prefix = "C:/Users/sheng/projects/paula_hurley/20201208_scRNA_split/split_samples/HYW_4701.HTO"
}else{
  h5file = args[1]
  output_prefix = args[2]
  hashtag_regex = args[3]
}

print(paste0("h5file=", h5file))
print(paste0("output_prefix=", output_prefix))
print(paste0("hashtag_regex=", hashtag_regex))

split(h5file, output_prefix, hashtag_regex)
