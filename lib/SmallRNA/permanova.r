
# permanova
library(vegan)
library(data.table)
library("RColorBrewer")
library(scales)

option_table<-read.table(parSampleFile4, sep="\t", stringsAsFactor=F)
myoptions<-split(option_table$V1, option_table$V2)
log_transform<-ifelse(myoptions$log_transform == "1", TRUE, FALSE)

if(is.null(myoptions$draw_nmds)){
  draw_nmds<-TRUE
}else{
  draw_nmds<-ifelse(myoptions$draw_nmds == "1", TRUE, FALSE)
}

file_list<-read.table(parSampleFile1, sep="\t", stringsAsFactors=F)
file_list$name<-tools::file_path_sans_ext(basename(file_list$V1))
file_list<-file_list[!grepl("_sequence$|minicontig$|isomiR$|read$|^tRNA_pm|^rRNA_pm|other$", file_list$name),]
file_list$name<-gsub("smallRNA_1mm_.+\\.","",file_list$name)
file_list$name<-gsub("bacteria_group1.+","Microbiome",file_list$name)
file_list$name<-gsub("bacteria_group2.+","Environment",file_list$name)
file_list$name<-gsub("_group.+","",file_list$name)

files <- file_list$V1

design_files<-read.table(parSampleFile3, sep="\t",  stringsAsFactors=F)

design_con <- lapply(design_files$V1, fread, data.table=F, select=c(1,2))
design_con <- do.call(rbind, design_con)
samples<- unique(design_con$Sample)
groups<-unique(design_con$Condition)
group_colors<-rainbow(length(groups))
names(group_colors)<-groups

librarySize<-read.csv(parFile3, row.names=1,check.names=FALSE,header=T,stringsAsFactors=F)
librarySize<-unlist(librarySize["Reads for Mapping",,drop=T])

###############################################################################################################################
summary_table_betadisper <- data.frame(matrix(NA, nrow = nrow(file_list), ncol = nrow(design_files)))
rownames(summary_table_betadisper)<-file_list$name
colnames(summary_table_betadisper)<-design_files$V2
summary_table_permanova <- summary_table_betadisper
summary_table_betadisper_f <- summary_table_betadisper
summary_table_permanova_f <- summary_table_betadisper

log_suffix<-ifelse(log_transform, ".log2", "")

if(draw_nmds){
  k_tbl <- NULL 
}

pdf(paste0(outFile, log_suffix, ".PCoA.pdf"))
i=1
for (i in 1:nrow(file_list)){
  count_name=file_list$name[i]
  count_file=file_list$V1[i]

  cat(count_name, "\n")

  data <- read.table(files[i], sep="\t", header=T, row.names=1, stringsAsFactors=F)
  data<-data[,colnames(data) %in% samples]

  #normalize total reads per million 
  count <- rbind(data, librarySize[colnames(data)])
  count_N <- data.frame(row.names=rownames(count),lapply(count, function(x)
    (x/x[nrow(count)])*1e+6))
  count_N <- count_N[-nrow(count_N),]
  
  count_N <- t(count_N)
  count_N <- count_N[sort(rownames(count_N)),]
  
  j=1
  for( j in 1:nrow(design_files)){
    design_name=design_files$V2[j]
    design_file=design_files$V1[j]
    cat("  ", design_name, "\n")

    prefix<-paste0(count_name, ".", design_name, log_suffix)

    design_con_sel<-read.table(design_files$V1[j], sep="\t", header=T, stringsAsFactors=F)
    design_con_sel<-design_con_sel[order(design_con_sel$Sample),]
    gcolors<-group_colors[unique(design_con_sel$Condition)]
    design_con_sel$color<-gcolors[design_con_sel$Condition]

    if (log_transform){
      count_N_sel<-log2(count_N[design_con_sel$Sample,]+1)
    }else{
      count_N_sel<-count_N[design_con_sel$Sample,]
    }

    dist <- vegdist(count_N_sel,  method = "bray")
    ordination <- wcmdscale(dist, eig = TRUE) 
    main_title = paste("PCoA", count_name, design_name, sep=" ")

    # Calculate variance explained
    Eigenvalues <- eigenvals(ordination) 
    Variance <- Eigenvalues / sum(Eigenvalues) 

    plot(ordination,type="n",  main = main_title, xlab=paste0("PCoA 1 (", label_percent(0.01)(Variance[1]),")"), ylab=paste0("PCoA 2 (", label_percent(0.01)(Variance[2]),")"))
    ordiellipse(ord=ordination, groups = design_con_sel$Condition, col=gcolors, kind="ehull", label=F)
    ordihull(ord=ordination, groups = design_con_sel$Condition, col=gcolors, display = "sites")
    ordispider(ord=ordination, groups = design_con_sel$Condition, col=gcolors, label =F)
    points(ordination$points[,1], ordination$points[,2], display = "sites", cex=1.5, pch=19,col=design_con_sel$color)
    legend("top", legend=unique(design_con_sel$Condition), col=gcolors, lty=c(1,1), bty = "n")

    if(draw_nmds){
      #https://cran.r-project.org/web/packages/vegan/vignettes/FAQ-vegan.html
      k <-1
      klim<-nrow(count_N_sel)-2
      bSucceed=FALSE
      error = ""
      while(!bSucceed & (k < klim)){
        ordination<-NA
        k = k+1
        cat(k, "\n")
        tryCatch(
          {
            ordination<-metaMDS(count_N_sel, distance="bray",trymax=20, autotransform = F, k=k, zerodist="add", weakties = FALSE)
            bSucceed=TRUE
          },
          warning=function(cond){
            error<<-paste0(error, "\nk=", k, ", warning: ", cond[[1]])
            cat("warning: ", cond[[1]], "\n")
            if(!is.na(ordination)){
              bSucceed<<-TRUE
            }
          },
          error=function(cond){
            error<<-paste0(error, "\nk=", k, ", error: ", cond[[1]])
            cat("error: ", cond[[1]], "\n")
          }
        )
      }

      main_title = paste("NMDS", count_name, design_name, sep=" ")

      if(is.na(ordination)){
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n', main = main_title)
        text(x = 0.5, y = 0.5, paste0(main_title, "\n\n",
                              "metaMDS failed due to error:\n",
                              error), 
            cex = 1, col = "black")      
      }else{
        plot(ordination,type="n",  main = main_title)
        ordihull(ord=ordination, groups = design_con_sel$Condition, col=gcolors, display = "sites")
        ordiellipse(ord=ordination, groups = design_con_sel$Condition, col=gcolors, kind="ehull", label=F)
        ordispider(ord=ordination, groups = design_con_sel$Condition, col=gcolors, label =F)
        points(ordination, display = "sites", cex=1.5, pch=19,col=design_con_sel$color)
        legend("top", legend=unique(design_con_sel$Condition), col=gcolors, lty=c(1,1), bty = "n")
      }

      if (is.null(k_tbl)){
        k_tbl<-data.frame("count"=count_name, "design"=design_name, "k"=k)
      }else{
        k_tbl<-rbind(k_tbl, data.frame("count"=count_name, "design"=design_name, "k"=k))
      }
    }

    #betadisper
    distance<-vegdist(count_N_sel)

    write.csv(as.matrix(distance), paste0(prefix, ".distance.csv"))

    betadisper <- anova(betadisper(distance, design_con_sel$Condition))
    summary_table_betadisper[i,j] <- betadisper[1,5]
    summary_table_betadisper_f[i,j] <- betadisper[1,4]

    #permanova
    permanova <- adonis(count_N_sel ~ design_con_sel[,2], method = "bray", perm =999)
    summary_table_permanova[i,j] <- permanova$aov.tab["Pr(>F)"][1,1]
    summary_table_permanova_f[i,j] <- permanova$aov.tab["F.Model"][1,1]
  }  
}
dev.off()

if(draw_nmds){
  write.csv(k_tbl, paste0(outFile, log_suffix, ".PCoA_K.csv"))
}

write.csv(summary_table_betadisper, paste0(outFile, log_suffix, ".betadisper.p.csv"))
write.csv(summary_table_betadisper_f, paste0(outFile, log_suffix, ".betadisper.f.csv"))
write.csv(summary_table_permanova, paste0(outFile, log_suffix, ".permanova.p.csv"))
write.csv(summary_table_permanova_f, paste0(outFile, log_suffix, ".permanova.f.csv"))
#end
