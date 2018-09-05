library(ggplot2)
library(data.table)

date <- gsub("-", "", Sys.Date())
#RNA_class <- c("miRNA","other","rRNA","snoRNA","snRNA","tRNA")

fileList1 <- fread("fileList1.txt", header=FALSE, data.table = F)
RNA_class <- fileList1[1:7,2]


category <- fread(fileList1[which(fileList1[,2]=="category"),1], data.table = F)
colnames(category)[1] <- "Sequence"



dirs <- fileList1[1:7, 1]

read.count.list<-list()
for (i in 1:length(dirs)){
  read.count<-fread(dirs[i],data.table=FALSE)
  colnames(read.count)[1]<-"Sequence"
  read.count<-rbind(read.count,category[5,])
  read.count<-data.frame(row.names=read.count[,1],lapply(read.count[2:ncol(read.count)],function(x)
    (x)/x[nrow(read.count)]*1e+6))
  read.count<-read.count[-nrow(read.count),]
  Sum_all<-data.frame(row.names=colnames(read.count[1:ncol(read.count)]))
  for (j in 16:67){
    Sum<-data.frame(colSums(read.count[which(nchar(rownames(read.count))==j),1:ncol(read.count)]))
    Sum_all<-cbind(Sum_all,Sum)}
  colnames(Sum_all)<-16:67
  read.count.list[[i]]<-Sum_all}


sample_list<-list()
for (i in 1:(ncol(category)-1)){
  sample_list[[i]]<-data.frame(rbind(read.count.list[[1]][i,],read.count.list[[2]][i,],
                                     read.count.list[[3]][i,],read.count.list[[4]][i,],
                                     read.count.list[[5]][i,],read.count.list[[6]][i,],
                                     read.count.list[[7]][i,]), 
                               row.names = RNA_class)
}

names(sample_list) <- rownames(read.count.list[[1]])

fastq_length<-fread(fileList1[9,1], data.table=FALSE)

fastq_length<-rbind(fastq_length[17:68,], category[5,2:ncol(category)])
fastq_length<-data.frame(row.names=16:68,lapply(fastq_length, function(x)
  (x/x[nrow(fastq_length)])*1e+6))
fastq_length<-fastq_length[-nrow(fastq_length),]
fastq_length<-t(fastq_length)


names(sample_list) <- rownames(read.count.list[[1]])

for (k in 1:length(sample_list)){
df<-data.frame(melt(sample_list[[k]]), RNA=factor(rep(c("miRNA","osDR","rDR","snoDR","snDR","tDR", "Genome"),52)
                          ,levels=rev(c("miRNA","osDR","rDR","snoDR","snDR","tDR", "Genome"))))
print(ggplot()+ geom_area(aes(x=16:67,y=fastq_length[k,]),fill="gray75",color="black")+
geom_bar(aes(x=rep(c(16:67),each=7), y=df[,2],fill=df[,3]),
         stat="identity", width=0.8, color="black",size=0.73)+
  scale_fill_manual(values=rev(c("blue","deeppink","red","purple","green","yellow", "black")))+
  theme(legend.title=element_blank())+
  theme_classic() + labs(x= NULL, y="Total Reads Per Million",title=names(sample_list)[k])+
  theme(plot.title = element_text(face= "bold", color = "black",size=22, hjust=0.5),axis.title = element_text(face = "bold", color = "black",size=20),
        axis.text = element_text(face= "bold", color = "black"),axis.text.y=element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.line = element_line(colour = "gray75", size =0.73, linetype = "solid"),
        axis.ticks = element_line(size=0.73),axis.ticks.length=unit(0.3,"cm"),
        legend.title=element_blank(), legend.position="none")+
  scale_y_continuous(expand = c(0, 0)) +scale_x_continuous(breaks=seq(16,67,by=4), labels=seq(16,67,by=4),expand = c(0, 0)))
ggsave(paste0(date, "_Lengh_BarPlot_Genome_", names(sample_list)[k], ".png"))
}
#END


