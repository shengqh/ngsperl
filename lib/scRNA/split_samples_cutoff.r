library(Seurat)
library(ggplot2)

#devtools::install_github("choisy/cutoff")
#install.packages("bbmle")
library(choisycutoff)
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
  den <- tryCatch(
    expr = {
      res=density(values, bw="SJ")
      return(res)
    },
    error = function(e){ 
      res=density(values)
      return(res)
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
  res=res[order(res$y, decreasing = T),]
  if(nrow(res)>2){
    #the highest peaks should be the negative one, positive one is always in the right side.
    res=res[res$x>=res$x[1],]

    if(nrow(res)>2){
      res=res[1:2,]
    }
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
      meta=t(meta)
    }
  }else{
    sdata<-Read10X_h5(h5file)
    meta<-sdata[[2]]
    meta<-as.matrix(meta)
  }
  mat<-meta
  write.csv(mat, file=paste0(output_prefix, ".alltags.exp.csv"))

  cat("All tags: ", paste(rownames(mat), collapse=","), "\n")
  
  if (!is.na(hashtag_regex)) {
    htos<-mat[grepl(hashtag_regex, rownames(mat)),]
    if (nrow(htos) == 0){
      stop(paste0("Cannot find hashtag based on regex ", hashtag_regex, " for tags ", paste(rownames(mat), collapse=",")))
    }
    cat("All hash tags matched regex: ", paste(rownames(htos), collapse=","), "\n")
  }else{
    htos<-mat
  }
  rownames(htos)<-gsub("^TotalSeqC_", "", rownames(htos))
  rownames(htos)<-gsub("^TotalSeq_", "", rownames(htos))
  rownames(htos)<-gsub('.TotalSeqC$', "", rownames(htos))

  empty_cell_sum<-apply(htos, 2, sum)
  htos<-htos[,empty_cell_sum > 0]

  write.csv(htos, file=paste0(output_prefix, ".hto.exp.csv"))
  
  pbmc.hashtag <- CreateSeuratObject(counts = htos, assay="HTO")
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")
  DefaultAssay(object = pbmc.hashtag) <- "HTO"
  
  #Idents(pbmc.hashtag) <- "HTO_classification"
  tagnames=rownames(pbmc.hashtag[["HTO"]])
  
  width=max(10, length(tagnames) * 5)
  pdf(paste0(output_prefix, ".tag.dist.pdf"), width=width, height=6)
  rplot(pbmc.hashtag, assay="HTO", features = tagnames, identName="orig.ident")
  dev.off()
  
  data <- FetchData(object=pbmc.hashtag, vars=tagnames)
  colnames(data)<-gsub("hto_","",colnames(data))

  tagname=tagnames[1]  
  for (tagname in tagnames) {
    values=data[,tagname]
    values=values[values>0]
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
  rplot(pbmc.hashtag, assay = "HTO", features = tagnames, identName="HTO_classification")
  dev.off()
  
  if (length(tagnames) == 2) {
    pdf(paste0(output_prefix, ".class.point.png"), width=width, height=height)
    print(FeatureScatter(object = pbmc.hashtag, feature1 = tagnames[1], feature2 = tagnames[2],group.by="HTO_classification"))
    dev.off()
  }
  
  tmat=data.frame(t(htos))
  tmat$HTO = obj$HTO_classification
  tmat$HTO.global = obj$HTO_classification.global
  write.csv(tmat, paste0(output_prefix, ".csv"))
  
  VariableFeatures(pbmc.hashtag)<-rownames(pbmc.hashtag)
  pbmc.hashtag<-ScaleData(pbmc.hashtag)
  pbmc.hashtag<-RunUMAP(pbmc.hashtag, features=rownames(pbmc.hashtag))
  
  png(paste0(output_prefix, ".umap.class.png"), width=1000, height=800)
  g<-DimPlot(pbmc.hashtag, reduction = "umap", group.by="HTO_classification")
  print(g)
  dev.off()
  
  png(paste0(output_prefix, ".umap.tag.png"), width=1600, height=1600)
  g<-FeaturePlot(pbmc.hashtag, features=tagnames, reduction = "umap")
  print(g)
  dev.off()
  
  hto_names=unique(pbmc.hashtag$HTO_classification)
  a_hto_names=hto_names[!(hto_names %in% c("Doublet","Negative"))]
  a_hto_names=a_hto_names[order(a_hto_names)]
  hto_names=c(a_hto_names, "Negative", "Doublet")
  cols=rep("gray", length(hto_names))
  names(cols)=hto_names
  cols[['Negative']]="blue"
  cols[["Doublet"]]="red"

  pbmc.hashtag$HTO_classification=factor(pbmc.hashtag$HTO_classification, levels=hto_names)
  png(paste0(output_prefix, ".umap.all.png"), width=1000, height=800)
  g<-DimPlot(pbmc.hashtag, reduction = "umap", label=T, group.by="HTO_classification", order=c("Negative", "Doublet"))+
    scale_color_manual(values=cols)
  print(g)
  dev.off()
}

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  h5file = r"(C:\projects\scratch\cqs\shengq2\paula_hurley_projects\20211130_scRNA_human\hto_samples_preparation\result\HYW_4701.hto.rds)"
  output_prefix = r"(C:\projects\scratch\cqs\shengq2\paula_hurley_projects\20211130_scRNA_human\HYW_4701.hto)"
  hashtag_regex='Hashtag|TotalSeqC_|C025|Benign|Tumor|HTO|HEK|THP|K562|KG1'
}else{
  h5file = args[1]
  output_prefix = args[2]
  hashtag_regex = args[3]
}

print(paste0("h5file=", h5file))
print(paste0("output_prefix=", output_prefix))
print(paste0("hashtag_regex=", hashtag_regex))

split(h5file, output_prefix, hashtag_regex)
