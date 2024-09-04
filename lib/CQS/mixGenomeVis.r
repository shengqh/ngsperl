rm(list=ls()) 
outFile='P11830_CW_combined'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/data/wanjalla_lab/projects/20240729_11830_RNAseq_mixedGenome/genetable/result/P11830_CW_combined.count'
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20240729_11830_RNAseq_mixedGenome/mix_genome_barplot/result')

### Parameter setting end ###

source("countTableVisFunctions.R")
library(data.table)
library(tibble)

options(bitmapType='cairo')

resultFile<-outFile

myoptions=read_file_map(parSampleFile1, do_unlist=FALSE)
host_prefix=gsub("_", "", myoptions$host_prefix)
cat("host_prefix=", host_prefix, "\n")

sample_groups=fread(parSampleFile2, header=FALSE)
sample_names=unique(sample_groups$V1)

counts = fread(parFile1, data.table=FALSE) |>
  dplyr::mutate(Genome=gsub("_.*","",Feature_chr)) |>
  dplyr::select(Genome, all_of(sample_names)) |>
  dplyr::group_by(Genome) |>
  dplyr::summarize(across(all_of(sample_names), sum))

gcounts = t(counts |> tibble::column_to_rownames("Genome")) |> as.data.frame()
gcounts$Total=rowSums(gcounts[,1:2])
gcounts$HostPercentage=gcounts[,host_prefix] / gcounts$Total * 100
gcounts$NonHostPercentage=100 - gcounts$HostPercentage
gcounts=gcounts |> tibble::rownames_to_column("Sample")

if(all(gcounts$Sample %in% sample_groups$V1)){
  gcounts$Group=sample_groups$V2[match(gcounts$Sample, sample_groups$V1)]
  gcounts=gcounts |> dplyr::arrange(Group, Sample)
  gcounts$Group=factor(gcounts$Group, levels=unique(gcounts$Group))
}  
gcounts$Sample=factor(gcounts$Sample, levels=gcounts$Sample)
write.csv(gcounts, paste0(resultFile,".Reads.csv"), row.names=FALSE)

if("Group" %in% colnames(gcounts)){
  g=ggplot(gcounts,aes(x=Sample,y=NonHostPercentage,fill=Group))
  height=5
}else{
  g=ggplot(gcounts,aes(x=Sample,y=NonHostPercentage))
  height=4
}
g=g+ylab("Non-host reads (%)") +
  geom_bar(stat="identity", width=0.5) +
  theme_bw2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        legend.position="top",
        legend.title=element_blank())

width<-min(max(14, length(unique(gcounts$Sample)) * 0.15), 60)
ggsave(paste0(outFile,".bar.png"),plot=g,width=width, height=height, dpi=300, units="in", bg="white")
ggsave(paste0(outFile,".bar.pdf"),plot=g,width=width, height=height)
