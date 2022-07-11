rm(list=ls()) 
outFile='AG3669'
parSampleFile1=''
parSampleFile2=''
parSampleFile3=''
parFile1='/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_4_consensus/result/metadata.txt'
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_5_immunarch/result')

### Parameter setting end ###

##Immunarch
library(immunarch) 
library(data.table)
library(knitr)

project_name <- outFile
target_dir <- paste0(getwd(), "/")

#Load
metadata <- parFile1
path <- dirname(metadata)
immdata <- repLoad(.path = path)
sample_names <- immdata$meta[,1]
plot_files <- c()

#Explore and compare T-cell and B-cell repertoires
repExplore(immdata$data, "lens") %>% vis()
cd3_dist <- paste0(target_dir, project_name, "_CD3_distribution.png")
ggsave(cd3_dist)

repClonality(immdata$data, "homeo") %>% vis() 
clon_summary <- paste0(target_dir, project_name, "_clonotype_summary.png")
ggsave(clon_summary)

repOverlap(immdata$data) %>% vis()
rep_overlap <- paste0(target_dir, project_name, "_repertoire_overlap.png")
ggsave(rep_overlap)

#Gene usage by sample
for (i in 1:nrow(immdata$meta)){
  sample_geneusage <- geneUsage(immdata$data[[i]]) 
  p <- sample_geneusage[which(sample_geneusage$Clones > 1), ] %>% vis() 
  p + ggtitle(paste0(sample_names[i,1], " Gene Usage"))
  p
  plot_name <- paste0(target_dir, sample_names[i,1], "_geneUsage.png")
  ggsave(plot_name)
  plot_files <- c(plot_files, plot_name)
}


repDiversity(immdata$data) %>% vis(.by = "Condition", .meta = immdata$meta)
diversity <- paste0(target_dir, project_name, "_diversity.png")
ggsave(diversity)

plot_files <- data.frame(c(plot_files, cd3_dist, clon_summary, rep_overlap, diversity, metadata))
write.table(plot_files, "plot_files.txt", col.names = F, row.names = F, quote=F, sep="\t")

output_file = paste0(project_name, ".html")
cat("Output report to:", target_dir, "\n")
rmarkdown::render("immunarch.Rmd",
                  output_dir = target_dir,
                  output_file = output_file)
##End
