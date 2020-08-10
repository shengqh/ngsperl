rm(list=ls())
library(readr)
library(stringr)
library(RColorBrewer)
library(pheatmap)
library(forcats)
library(dendextend)
library(getopt)
library(Cairo)
spec =matrix(c(
  'output','o', 1,"character",
  'inputscore','i', 1,"character",
  'inputfamily','f', 1,"character"
), byrow=TRUE, ncol=4)
opt =getopt(spec)
output<-opt$output
inputscore<-opt$inputscore
inputfamily<-opt$inputfamily
my_file<-inputscore
family_info_file<-inputfamily

#myfile<-"C:/Users/lenovo/Desktop/ibs/Ciombor_ExomeSeq.ibs_score.mean.csv"
#family_info_file<-"C:/Users/lenovo/Desktop/ibs/Ciombor_ExomeSeq__fileList1.list"
method<-"mean"
data<-read.csv(my_file, header = TRUE, sep = ",",row.names = 1,check.names = F)

family_info<-read.table(family_info_file,sep = "\t",header = F,check.names = F,colClasses = c("character","factor"))
colnames(family_info)<-c("sampleID","familyID")

data[lower.tri(data,diag=T)]<-t(data)[lower.tri(data,diag=T)]

### add familyID to sampleID in data
data_sname<-data.frame(sampleID=colnames(data),stringsAsFactors = F)
sample_finfo<-merge(data_sname,family_info,by = "sampleID",sort = F)
data$familyID<-sample_finfo$familyID
########### order data in the order of familyID factor
data<-data[order(data$familyID),]
data<-as.data.frame(t(data[,-1*ncol(data)]),stringsAsFactors = F)
data$familyID<-sample_finfo$familyID
data<-data[order(data$familyID),]
data<-data[,-1*ncol(data)]
####### delete the rows and cols with fulling with only NA
if(length(which(apply(data,1,function(x)all(is.na(x)))))>0)
{
  data<-data[-which(apply(data,1,function(x)all(is.na(x)))), -which(apply(data,2,function(x)all(is.na(x))))]
}
#### set the color in plot
cc = colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))

annotation_col<-data.frame(familyID=factor(sample_finfo$familyID[sapply(colnames(data),function(x) which(sample_finfo$sampleID==x))]),row.names = colnames(data))
#ann_colors = list( GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E"))
#heatmap=pheatmap(data,color = cc(100),
#                main="IBS",
#                fontsize = 1000/ncol(data),
#                 #scale="row",
#                 border_color = NA,
#                 na_col = "grey",
#                 cluster_rows = F,cluster_cols = F,
#                 show_rownames = T,show_colnames = T,
#                treeheight_row = 30,treeheight_col = 30,
#                 cellheight = 1000/ncol(data),cellwidth = 1000/ncol(data),
#                 #cutree_row=2,cutree_col=2,
#                 display_numbers = F,legend = T,
#                 annotation_col = annotation_col,
#                 #annotation_colors = ann_colors,
#                 #annotation_legend = FALSE,
#                 filename = paste0("IBS_",method,"_uncluster.png")
#)

dend <- as.dendrogram(hclust(2-as.dist(data)))
cluster_sample<-data.frame(sampleID=labels(dend),stringsAsFactors = F)

#### record cluster result with the frequency of familyID (the number of members) in the cluster
annotation_col$sampleID<-rownames(annotation_col)
cluster_result<-merge(cluster_sample,annotation_col,by="sampleID",sort=F)
countmembers<-data.frame(fct_count(cluster_result$familyID))
cluster_result$n_members<-countmembers$n[sapply(as.character(cluster_result$familyID),function(x) which(as.character(countmembers$f)==x))]
cluster_result$familyID<-as.character(cluster_result$familyID)
### record the first and last index of sampleID in the cluster to check whether cluster together or not
check_keep_dup<-data.frame(familyID=unique(cluster_result$familyID),stringsAsFactors = F)
sindex<-sapply(check_keep_dup$familyID,function(x) which(cluster_result$familyID==x),simplify = F)
for(i in 1:NROW(sindex)) {check_keep_dup$first_index[i]=sindex[[i]][1]}
for(i in 1:NROW(sindex)) {check_keep_dup$last_index[i]=sindex[[i]][length(sindex[[i]])]}
###  check whether cluster together or not

cluster_result<-merge(cluster_result,check_keep_dup,by = "familyID",sort = T)
# order the cluster_result with the samppleID in the order of cluster
cluster_result<-merge(cluster_sample,cluster_result,by="sampleID",sort = F)
cluster_result$cluster_together<-(cluster_result$first_index+cluster_result$n_members-1)==cluster_result$last_index
cluster_result$cluster_together<-ifelse(cluster_result$cluster_together,"YES","NO")
cluster_result$part_cluster_together<-NA
ind<-1
while(ind<nrow(cluster_result))
{
  if(cluster_result$familyID[ind]==cluster_result$familyID[ind+1])
  {
    cluster_result$part_cluster_together[ind]<-"YES"
    cluster_result$part_cluster_together[ind+1]<-"YES"
    #ind<-ind+1
  }
  else
  {
    if(is.na(cluster_result$part_cluster_together[ind]))
      cluster_result$part_cluster_together[ind]<-"NO"
    if(!is.na(cluster_result$part_cluster_together[ind]))
    {
      if(cluster_result$part_cluster_together[ind]=="YES")
        cluster_result$part_cluster_together[ind]<-"YES"
      else
      {
        cluster_result$part_cluster_together[ind]<-"NO"
      }
    }
    cluster_result$part_cluster_together[ind+1]<-"NO"
  }
  ind<-ind+1
}
cluster_result$part_cluster_together[which(cluster_result$cluster_together=="YES")]<-"YES"

cluster<-cluster_result$cluster_together


cluster[which(cluster_result$cluster_together=="YES"&cluster_result$part_cluster_together=="YES")]<-"red"

cluster[which(cluster_result$cluster_together=="NO"&cluster_result$part_cluster_together=="YES")]<-"darkorange2"

cluster[which(cluster_result$cluster_together=="NO"&cluster_result$part_cluster_together=="NO")]<-"moccasin"
cluster_result$cluster<-cluster

the_bars <- data.frame(cluster=cluster)



CairoPNG(paste0(output, "_cluster.png"),width=ncol(data)*25,height=ncol(data)*25,res=ncol(data)*3)
par(mar = c(2,1,4,12))
dend %>% set("labels_cex", 1000/ncol(data)/30) %>% 
plot(main = "IBS",horiz = T) # change color 
colored_bars(colors = the_bars, dend = dend, sort_by_labels_order = F,horiz = T, y_shift= 1000/ncol(data)/30,rowLabels="")

#colored_bars(colors = the_bars, dend = dend, sort_by_labels_order = F,horiz = T,y_shift= 1000/ncol(data)/30*5, rowLabels="")

legend("topleft", legend = c("Whole_Family_Members", "Part_Family_Members","Seperate_Family_Members"), pch = 15, cex = 1000/ncol(data)/30*5, bty = 'n',
       #inset = c(1, 1), # place outside
       title = "cluster", 
       col = c("red","darkorange2","moccasin"))

dev.off()
write.table(cluster_result[,c(1:3,6:8)],paste0(output, "_cluster.csv"),sep = ",",row.names =F, quote =F)

