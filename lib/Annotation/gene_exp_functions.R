#### Functions for my section ####
calc_z<-function(x, u, s){
  z = (x-u)/s
  return(z)
}

tiss_GT1_z<-function(to_plot_zScored){
  list_of_tiss_genes<-c()
  for(i in 1:ncol(to_plot_zScored)){
    list_of_tiss_genes<-c(list_of_tiss_genes,list(row.names(to_plot_zScored)[which(to_plot_zScored[,i]>1)]))
  }
  names(list_of_tiss_genes)<-names(to_plot_zScored)
  
  counts_of_tiss_genes<-as.data.frame(matrix(nrow = length(list_of_tiss_genes), ncol = 2))
  names(counts_of_tiss_genes)<-c("Tissue","NumGenesGT1_zscore")
  counts_of_tiss_genes$Tissue<-names(list_of_tiss_genes)
  counts_of_tiss_genes$NumGenesGT1_zscore<-unlist(lapply(list_of_tiss_genes, length))
  
  counts_ordered<-counts_of_tiss_genes[order(counts_of_tiss_genes$NumGenesGT1_zscore, decreasing = F),]
  
  counts_ordered$Tissue<-factor(counts_ordered$Tissue, levels = unique(counts_ordered$Tissue))
  
  return(counts_ordered)
}

proteins_GT1_zscore_diff_tissNum<-function(to_plot_zScored, proteins){
  protein_counts<-c()
  for(i in 1:nrow(to_plot_zScored)){
    protein_counts<-c(protein_counts, length(which(to_plot_zScored[i,]>1)))
  }
  names(protein_counts)<-row.names(to_plot_zScored)
  
  protein_counts_plot<-as.data.frame(matrix(nrow = length(protein_counts), ncol = 3))
  names(protein_counts_plot)<-c("Ensembl","Gene_Sym","counts")
  protein_counts_plot$Ensembl<-row.names(to_plot_zScored)
  protein_counts_plot$Gene_Sym<-proteins$GENE[match(row.names(to_plot_zScored), proteins$ensemblID)]
  protein_counts_plot$counts<-unlist(protein_counts)
  
  return(protein_counts_plot)
}

top_5perc_genes<-function(red_tabula){
  top5perc<-list()
  for(i in 1:ncol(red_tabula)){
    
    the_geneExpr<-red_tabula[which(red_tabula[,i]>0),i]
    the_geneExpr_genes<-row.names(red_tabula)[which(red_tabula[,i]>0)]
    
    the_ranks<-order(the_geneExpr, decreasing = T)
    the_gene_set<-the_geneExpr_genes[which(the_ranks<(0.05*length(the_geneExpr)))]
    
    top5perc[[i]]<-list(the_gene_set)
  }
  
  names(top5perc)<-names(red_tabula)
  
  return(top5perc)
}

counts_and_ratio_top5perc<-function(top5perc, red_proteins, pval, log_for_filt, target_column, GeneID1, GeneID2){
  cou_rat_tab<-as.data.frame(matrix(nrow = length(top5perc), ncol = 6))
  names(cou_rat_tab)<-c("Tabula_group","Count_sig_proteins_in_top5perc","Num_of_top5perc_genes","Ratio_sigProteins_in_NumTop5PercGenes", "Proteins_geneID1","Proteins_geneID2")
  cou_rat_tab$Tabula_group<-names(top5perc)
  
  sig_proteins<-red_proteins[which(red_proteins[,target_column]),]
  
  for(i in 1:nrow(cou_rat_tab)){
    if(table(sig_proteins$ensemblID%in%unlist(top5perc[[i]]))['FALSE']!=nrow(sig_proteins)){
      cou_rat_tab$Count_sig_proteins_in_top5perc[i]<-table(unlist(top5perc[[i]])%in%sig_proteins[,GeneID2])['TRUE']
      cou_rat_tab$Proteins_geneID1[i]<-paste(sig_proteins[which(sig_proteins$ensemblID %in% unlist(top5perc[[i]])), GeneID1], sep = ",", collapse = ",")
      cou_rat_tab$Proteins_geneID2[i]<-paste(sig_proteins[which(sig_proteins$ensemblID %in% unlist(top5perc[[i]])), GeneID2], sep = ",", collapse = ",")
    }else{
      cou_rat_tab$Count_sig_proteins_in_top5perc[i]<-0
    }
    
    cou_rat_tab$Num_of_top5perc_genes[i]<-length(unlist(top5perc[[i]]))
    
    cou_rat_tab$Ratio_sigProteins_in_NumTop5PercGenes[i]<-(cou_rat_tab$Count_sig_proteins_in_top5perc[i]/cou_rat_tab$Num_of_top5perc_genes[i])
  }
  
  return(cou_rat_tab)
}

cumulativePlot=function(count_matrix, order_method = "max", sample_to_use = NULL, label_to_use = NULL) {
  
  #order by
  if(length(order_method)==nrow(count_matrix)){
    orderBy = order_method
    if(is.null(label_to_use)){
      order_method = "vector_supplied_by_user"
    }else{
      order_method = label_to_use
    }
  }else{
    if(order_method=="max"){
      orderBy<-apply(count_matrix, 1, max)
    }
    if(order_method=="mean"){
      orderBy<-apply(count_matrix, 1, mean)
    }
    if(order_method=="median"){
      orderBy<-apply(count_matrix, 1, median)
    }
    if(order_method=="sample"){
      if(is.null(sample_to_use)){
        stop("Please supply the sample name to order by if you want to use sample ordering")
      }
      orderBy<-count_matrix[,sample_to_use]
      order_method = sample_to_use
    }
  }
  
  geneOrder=rev(row.names(count_matrix)[order(orderBy)])
  count_matrix=count_matrix[geneOrder,]
  
  cumulative_sum <-(apply(count_matrix, 2, cumsum))
  cumulative_percent <- apply(cumulative_sum, 2, function(x) x/max(x) * 100)
  #make sure things looks right
  head(cumulative_percent)
  
  # Convert matrix to data frame
  df <- as.data.frame(cumulative_percent)
  
  # Reshape data frame for ggplot2
  df$Gene <- row.names(df)
  df_long <- df %>% gather(Sample, Percentage, -Gene)
  
  df_long$Gene=factor(df_long$Gene,levels = geneOrder)
  
  #Label the top two samples
  min_cum_in_top_genes<-c()
  for(tiss in unique(df_long$Sample)){
    min_cum_in_top_genes<-c(min_cum_in_top_genes, min(df_long$Percentage[df_long$Sample==tiss][-c(1:20)]))
  }
  names(min_cum_in_top_genes)<-unique(df_long$Sample)
  
  which_label = names(min_cum_in_top_genes[order(min_cum_in_top_genes, decreasing = T)])[1:3]
  
  # Plot
  p=ggplot(df_long, aes(x = Gene, y = Percentage, group =Sample , color = Sample)) +
    geom_line() +
    labs(x = paste0("Gene ranked by ", order_method), y = "Cumulative Percentage", title = "Cumulative Percentage Plot") +
    theme_minimal() +
    annotate('text', x = 20, y = df_long$Percentage[which(df_long$Sample==which_label[1])][20], label = as.character(which_label[1]))+
    annotate('text', x = 25, y = df_long$Percentage[which(df_long$Sample==which_label[2])][25], label = as.character(which_label[2]))+
    annotate('text', x = 30, y = df_long$Percentage[which(df_long$Sample==which_label[3])][30], label = as.character(which_label[3]))
  # theme(legend.position = "none")
  return(p)
}



draw_barplot<-function(rankData, genes, fileprefix, width, height, score_name, top=50){
  scoredf <- simpleScore(rankData, upSet = genes)
  scoredf <- scoredf[order(scoredf$TotalScore, decreasing=T),]
  scoredf <- scoredf[1:min(nrow(scoredf), top),]
  scoredf$Tissue=factor(rownames(scoredf), levels=rownames(scoredf))
  
  
  g<-ggplot(scoredf, aes(Tissue, TotalScore)) + geom_bar(stat="identity") + theme_bw3() + xlab("") + ylab(score_name)
  
  pdf(paste0(fileprefix, ".bar.pdf"), width=width, height=height)
  print(g)
  dev.off()
  print(g)
}

draw_alt_barplot<-function(scoredf, fileprefix, width, height, score_name){
  g<-ggplot(scoredf, aes(Tissue, TotalScore)) + geom_bar(stat="identity") + theme_bw3() + xlab("") + ylab(score_name)
  
  pdf(paste0(fileprefix, ".bar.pdf"), width=width, height=height)
  print(g)
  dev.off()
  
  print(g)
}

cumulativePlot=function(count_matrix, order_method = "max", sample_to_use = NULL, label_to_use = NULL) {
  
  #order by
  if(length(order_method)==nrow(count_matrix)){
    orderBy = as.numeric(unlist(order_method))
    if(is.null(label_to_use)){
      order_method = "vector_supplied_by_user"
    }else{
      order_method = label_to_use
    }
  }else{
    if(order_method=="max"){
      orderBy<-apply(count_matrix, 1, max)
    }
    if(order_method=="mean"){
      orderBy<-apply(count_matrix, 1, mean)
    }
    if(order_method=="median"){
      orderBy<-apply(count_matrix, 1, median)
    }
    if(order_method=="sample"){
      if(is.null(sample_to_use)){
        stop("Please supply the sample name to order by if you want to use sample ordering")
      }
      orderBy<-count_matrix[,sample_to_use]
      order_method = sample_to_use
    }
  }
  
  geneOrder=rev(row.names(count_matrix)[order(orderBy)])
  count_matrix=count_matrix[geneOrder,]
  
  cumulative_sum <-(apply(count_matrix, 2, cumsum))
  cumulative_percent <- apply(cumulative_sum, 2, function(x) x/max(x) * 100)
  
  # Convert matrix to data frame
  df <- as.data.frame(cumulative_percent)
  
  # Reshape data frame for ggplot2
  df$Gene <- row.names(df)
  df_long <- df %>% gather(Sample, Percentage, -Gene)
  
  df_long$Gene=factor(df_long$Gene,levels = geneOrder)
  
  #Label the top two samples
  min_cum_in_top_genes<-c()
  for(tiss in unique(df_long$Sample)){
    min_cum_in_top_genes<-c(min_cum_in_top_genes, min(df_long$Percentage[df_long$Sample==tiss][-c(1:20)]))
  }
  names(min_cum_in_top_genes)<-unique(df_long$Sample)
  
  which_label = names(min_cum_in_top_genes[order(min_cum_in_top_genes, decreasing = T)])[1:3]
  
  # Plot
  p=ggplot(df_long, aes(x = Gene, y = Percentage, group =Sample , color = Sample)) +
    geom_line() +
    labs(x = paste0("Gene ranked by ", order_method), y = "Cumulative Percentage", title = "Cumulative Percentage Plot") +
    theme_minimal() +
    annotate('text', x = 20, y = df_long$Percentage[which(df_long$Sample==which_label[1])][20], label = as.character(which_label[1]))+
    annotate('text', x = 25, y = df_long$Percentage[which(df_long$Sample==which_label[2])][25], label = as.character(which_label[2]))+
    annotate('text', x = 30, y = df_long$Percentage[which(df_long$Sample==which_label[3])][30], label = as.character(which_label[3]))
  # theme(legend.position = "none")
  return(p)
}
