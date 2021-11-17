
# permanova
library(vegan)
library(data.table)
library("RColorBrewer")

file_list<-read.table(parSampleFile1, sep="\t")
file_list$name<-tools::file_path_sans_ext(basename(file_list$V1))
file_list<-file_list[!grepl("_sequence$|minicontig$|isomiR$|read$|^tRNA_pm|^rRNA_pm|other$", file_list$name),]
file_list$name<-gsub("smallRNA_1mm_.+\\.","",file_list$name)
file_list$name<-gsub("bacteria_group1.+","Microbiome",file_list$name)
file_list$name<-gsub("bacteria_group2.+","Environment",file_list$name)
file_list$name<-gsub("_group.+","",file_list$name)

files <- file_list$V1

design_files<-read.table(parSampleFile3, sep="\t")

design_con <- lapply(design_files$V1, fread, data.table=F, select=c(1,2))
design_con <- do.call(rbind, design_con)
samples<- unique(design_con$Sample)
groups<-unique(design_con$Condition)
group_colors<-rainbow(length(groups))
names(group_colors)<-groups

librarySize<-read.csv(parFile3, row.names=1,check.names=FALSE,header=T,stringsAsFactor=F)
librarySize<-unlist(librarySize["Reads for Mapping",,drop=T])

###############################################################################################################################
summary_table_betadisper <- data.frame(matrix(NA, nrow = nrow(file_list), ncol = nrow(design_files)))
rownames(summary_table_betadisper)<-file_list$name
colnames(summary_table_betadisper)<-design_files$V2
summary_table_permanova <- summary_table_betadisper

pdf(paste0(outFile, ".PCoA.pdf"))
i=6
for (i in 1:nrow(file_list)){
  count_name=file_list$name[i]
  count_file=file_list$V1[i]

  cat(count_name, "\n")

  data <- read.table(files[i], sep="\t", header=T, row.names=1)
  data<-data[,colnames(data) %in% samples]

  #normalize total reads per million 
  count <- rbind(data, librarySize[colnames(data)])
  count_N <- data.frame(row.names=rownames(count),lapply(count, function(x)
    (x/x[nrow(count)])*1e+6))
  count_N <- count_N[-nrow(count_N),]
  
  count_N <- t(count_N)
  count_N <- count_N[sort(rownames(count_N)),]
  
  j=4
  for( j in 1:nrow(design_files)){
    design_name=design_files$V2[j]
    design_file=design_files$V1[j]
    cat("  ", design_name, "\n")

    prefix<-paste0(count_name, ".", design_name)

    design_con_sel<-read.table(design_files$V1[j], sep="\t", header=T)
    design_con_sel<-design_con_sel[order(design_con_sel$Sample),]
    gcolors<-group_colors[unique(design_con_sel$Condition)]
    design_con_sel$color<-gcolors[design_con_sel$Condition]

    count_N_sel<-count_N[design_con_sel$Sample,]

    #https://cran.r-project.org/web/packages/vegan/vignettes/FAQ-vegan.html
    k <-max(2, floor((nrow(count_N_sel) - 1)/4))
    original_k=k
    bSucceed=FALSE
    while(!bSucceed){
      cat(k, "\n")
      tryCatch(
        {
          ordination<-metaMDS(count_N_sel, distance="bray",trymax=20, autotransform = T, k=k)
          bSucceed=TRUE
        },
        warning=function(cond){
          cat("warning: ", cond[[1]], "\n")
          k=k-1
          if(k < 2){
            stop('cannot find K for metaMDS')
          }
        },
        error=function(cond){
          cat("error: ", cond[[1]], "\n")
          k=k-1
          if(k < 2){
            stop('cannot find K for metaMDS')
          }
        }
      )
    }

    if(original_k != k){
      main_title = paste("PCoA", count_name, design_name, paste0("k=",k), sep=" ")
    }else{
      main_title = paste("PCoA", count_name, design_name, sep=" ")
    }
    plot(ordination,type="n",  main = main_title)
    ordihull(ord=ordination, groups = design_con_sel$Condition, col=gcolors, display = "sites")
    ordiellipse(ord=ordination, groups = design_con_sel$Condition, col=gcolors, kind="ehull", label=F)
    ordispider(ord=ordination, groups = design_con_sel$Condition, col=gcolors, label =F)
    points(ordination, display = "sites", cex=1.5, pch=19,col=design_con_sel$color)
    #text(ordination, display = "sites", cex=0.7)

    #legend
    legend("top", legend=unique(design_con_sel$Condition), col=gcolors, lty=c(1,1), bty = "n")

    #betadisper
    distance<-vegdist(count_N_sel)

    write.csv(as.matrix(distance), paste0(prefix, ".distance.csv"))

    betadisper <- anova(betadisper(distance, design_con_sel$Condition))
    summary_table_betadisper[i,j] <- betadisper[1,5]

    #permanova
    permanova <- adonis(count_N_sel ~ design_con_sel[,2], method = "bray", perm =999)
    summary_table_permanova[i,j] <- permanova$aov.tab["Pr(>F)"][1,1]
  }  
}
dev.off()

write.csv(summary_table_betadisper, paste0(outFile, ".betadisper.csv"))
write.csv(summary_table_permanova, paste0(outFile, ".permanova.csv"))
#end
