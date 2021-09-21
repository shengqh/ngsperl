
library(Seurat)
library(patchwork)
library(ggplot2)

souporcell_tb<-read.table(parSampleFile1, sep="\t", row.names=2)
cutoff_tb<-read.table(parSampleFile2, sep="\t", row.names=2)
umap_rds_tb<-read.table(parSampleFile3, sep="\t", row.names=2)

sample_name=rownames(souporcell_tb)[2]
for (sample_name in rownames(souporcell_tb)){
  s1<-read.csv(cutoff_tb[sample_name, "V1"], row.names=1)
  s2<-read.table(souporcell_tb[sample_name, "V1"], row.names=1, header=T)
  obj<-readRDS(umap_rds_tb[sample_name, "V1"])
  
  s <- merge(s1, s2, by=0, all.x=TRUE)
  
  single<-s[s$HTO.global == "Singlet",]
  single<-single[single$status == "singlet",]
  tb<-table(single$HTO, single$assignment)
  
  write.csv(tb, paste0(sample_name, ".HTO_soupor_singlet.csv"))
  
  maxrow<-unlist(rownames(tb)[apply(tb,2,which.max)])
  
  cmap<-split(maxrow, colnames(tb))

  final<-unlist(apply(s, 1, function(x){
    #print(x['HTO.global'])
    if (x['HTO.global'] == 'Singlet'){
      return(x['HTO'])
    }
    if(x['HTO.global'] == 'Negative'){
      return(x['HTO'])
    }
    if(x['status'] == 'doublet'){
      return(x['HTO'])
    }
    if(x['status'] == 'unassigned'){
      return(x['HTO'])
    }
    
    ass=x['assignment']
    
    if(!(ass %in% names(cmap))){
      return(x['HTO'])
    }
    
    res=cmap[ass]
    #print(paste0(x['assignment'], " => ", res))
    return(res)
  }))
  
  s$souporcell_cutoff=unlist(apply(s, 1, function(x){
    if(x['status'] == 'doublet'){
      return("Doublet")
    }
    if(x['status'] == 'unassigned'){
      return("Unassigned")
    }
    
    ass=x['assignment']
    
    if(!(ass %in% names(cmap))){
      return('unmapped')
    }
    
    res=cmap[ass]
    #print(paste0(x['assignment'], " => ", res))
    return(res)
  }))
  
  s$souporcell=unlist(apply(s, 1, function(x){
    if(x['status'] == 'doublet'){
      return("Doublet")
    }
    if(x['status'] == 'unassigned'){
      return("Unassigned")
    }
    
    ass=x['assignment']
    return(ass)
  }))
  
  s$Final<-final
  rownames(s)<-s$Row.names

  obj<-subset(obj, cells=rownames(s))
  s<-s[colnames(obj),]
  rownames(s)<-s$Row.names

  tags<-unlist(cmap)
  tags<-tags[order(tags)]
  tags<-gsub("-",".",tags)
  hto_final<-s[,c(tags, "Final")]
  hto_final$HTO.golbal=unlist(lapply(hto_final$Final, function(x){
    if(x == "Doublet"){
      return("Doublet")
    }
    if(x == "Negative"){
      return("Negative")
    }
    return("Singlet")
  }))
  colnames(hto_final)<-c(tags, "HTO", "HTO.golbal")
  write.csv(hto_final, paste0(sample_name, ".HTO.csv"))
  
  obj$souporcell<-s$souporcell
  obj$souporcell_cutoff<-s$souporcell_cutoff
  obj$final<-s$Final
  obj$assignment<-s$assignment

  pt.size=0.3

  ss<-s[s$status=="singlet",]
  ss_assign<-unique(ss$assignment)
  ss_assign<-ss_assign[order(ss_assign)]
  for (assign in ss_assign){
    g2<-DimPlot(obj, group.by="assignment", pt.size = pt.size) + ggtitle(paste0("soupor cluster ", assign))
    gdata<-g2$data
    gdata$assignment<-as.character(gdata$assignment)
    gdata$assignment[gdata$assignment != assign] <- "Other"
    gdata$assignment[gdata$assignment == assign] <- unlist(s[rownames(gdata)[gdata$assignment==assign], "HTO"])
    
    groups<-unique(as.character(gdata$assignment))
    groups<-groups[order(groups)]
    col1<-c(rainbow(length(groups)-1), "grey")
    names(col1)<-c(groups[groups != "Other"], "Other")
    g2<-g2+scale_color_manual(values=col1)
    g2$data<-rbind(gdata[gdata$assignment == "Other",], gdata[gdata$assignment != "Other",])
    g2$data[,3]=factor(g2$data[,3], levels=groups)

    png(paste0(sample_name, ".soupor_cluster", assign, ".png"), width=1500, height=1000, res=300)
    print(g2)
    dev.off()
  }
  
  g1<-DimPlot(obj, group.by="HTO_classification", pt.size = pt.size)
  g2<-DimPlot(obj, group.by="souporcell", pt.size = pt.size)
  g3<-DimPlot(obj, group.by="final", pt.size = pt.size)

  g1levels<-unique(as.character(g1$data[,3]))
  g1levels<-g1levels[order(g1levels)]
  g3levels<-unique(as.character(g3$data[,3]))
  g3levels<-g3levels[order(g3levels)]
  
  glevels=c(g1levels, g3levels)
  glevels<-unique(glevels)
  glevels<-glevels[order(glevels)]
  gcolors<-rainbow(length(glevels))
  names(gcolors)<-glevels
    
  g1colors<-gcolors[g1levels]
  g3colors<-gcolors[g3levels]
  
  g1<-g1+scale_color_manual(values=g1colors)
  g3<-g3+scale_color_manual(values=g3colors)
  
  g<-g1+g2+g3
  png(paste0(sample_name, ".HTO.png"), width=4000, height=1000, res=300)
  print(g)
  dev.off()
  
  get_gg<-function(g1, hto){
    groups<-unique(as.character(g1$data[,3]))
    groups<-groups[order(groups)]
    col1<-rep("gray", length(groups))
    names(col1)<-groups
    col1[hto]<-"red"
    gg1<-g1+scale_color_manual(values=col1)
    gg1$data<-rbind(gg1$data[gg1$data[,3] != hto,], gg1$data[gg1$data[,3] == hto,])
    gg1$data[,3]=factor(gg1$data[,3], levels=groups)
    return(gg1)
  }

  hmap<-split(names(cmap), as.character(cmap))  
  hmap$Doublet<-"Doublet"

  htos<-unique(as.character(obj$HTO_classification))
  hto<-htos[1]
  for (hto in htos){
    gg1<-get_gg(g1, hto)
    
    if(hto %in% names(hmap)){
      gg2<-get_gg(g2, as.character(hmap[hto]))
    }else{
      gg2<-get_gg(g2, "Unmapped")
    }
    gg3<-get_gg(g3, hto)
    gg<-gg1+gg2+gg3
    
    png(paste0(sample_name, ".", hto, ".png"), width=4000, height=1000, res=300)
    print(gg)
    dev.off()
  }
}
