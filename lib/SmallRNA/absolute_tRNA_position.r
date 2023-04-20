rm(list=ls()) 
outFile='9074_ES'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/vickers_lab/projects/20221122_9074_ES_ARMseq_human_byMars_tRNA_pos/host_tRNA_abs_position_vis/result')

### Parameter setting end ###

read_file_map<-function(file_list_path, sep="\t", do_unlist=TRUE, reverse=FALSE){
  tbl<-read.table(file_list_path, sep=sep, header=F)
  if(reverse){
    result<-split(tbl$V2, tbl$V1)
  }else{
    result<-split(tbl$V1, tbl$V2)
  }
  if(do_unlist){
    result<-unlist(result)
  }
  return(result)
}


theme_rotate_x_axis_label <- function(angle=90, vjust=0.5, hjust=1) {
  theme(axis.text.x = element_text(angle=angle, vjust=vjust, hjust=hjust))
}

theme_bw3 <- function (axis.x.rotate=F, angle=90, vjust=0.5, hjust=1) { 
  is_ggplot2_newver = packageVersion("ggplot2") >= "3.4.0"

  if(is_ggplot2_newver){
    eline = element_line(colour = "black", linewidth = 0.5)
  }else{
    eline = element_line(colour = "black", size = 0.5)
  }

	result = theme_bw() +
    theme(
      strip.background = element_rect(fill = NA, colour = 'black'),
      panel.border = element_rect(fill = NA, color = "black"),			
      axis.line = eline
    )
  if (axis.x.rotate){
    result = result + theme_rotate_x_axis_label(angle=angle,vjust=vjust,hjust=hjust)
  }

  return(result)
}


myoptions = read_file_map(parSampleFile4, do_unlist=FALSE)

position_map = read_file_map(parSampleFile1)

group_map = read_file_map(parSampleFile2, reverse=TRUE)

if(myoptions$sample_pattern != ""){
  group_map=group_map[grepl(myoptions$sample_pattern, names(group_map))]
}
position_map=position_map[names(group_map)]

all_positions = NULL
sample=names(position_map)[1]
for(sample in names(position_map)){
  cat("Reading ", sample, "...\n")
  filepath=position_map[sample]
  cur_positions = read.table(filepath, sep="\t", header=T)
  cur_positions$Sample = sample
  all_positions = rbind(all_positions, cur_positions)
}

all_positions$Category=factor(all_positions$Category, levels=c("Not Amino Specific", "Isotype Specific", "Isodecoder Specific", "Transcript Specific"))

all_positions$Group = group_map[all_positions$Sample]

library(dplyr)
group_positions = all_positions %>% group_by(Group, Transcript, Category, Position) %>% summarise(SumCount=sum(Count))
group_positions$Isotype = gsub("-[^-]+-[^-]+$", "", group_positions$Transcript)

library(ggplot2)
isotype = group_positions$Isotype[1]
for(isotype in unique(group_positions$Isotype)){
  cat("Processing", isotype, "...\n")
  isotype_positions = subset(group_positions, Isotype==isotype)

  g<-ggplot(isotype_positions, aes(Position, SumCount, fill=Category)) + 
    geom_bar(stat="Identity") + 
    facet_grid(Transcript~Group, scale="free_y") + 
    theme_bw3() +
    theme(strip.text.y.right = element_text(angle = 0, hjust = 0))

  n_transcript = length(unique(isotype_positions$Transcript))
  n_group = length(unique(isotype_positions$Group))

  width=max(2000, n_group * 600 + 1000)
  height=max(2000, n_transcript * 200)

  png(paste0(isotype, ".tRNA.abs_position.png"), width=width, height=height, res=300)
  print(g)
  dev.off()
}
cat("Done\n")
