dmcpgs <- function(file,
                   sample = NULL,
                   perc_cut = 0.25,
                   FDR = 0.05,
                   mincov = 4){
  lst <- list()
  samples <- unlist(strsplit(sample, ","))
  df <- read.table(file, header = F, sep = "\t")
  colnames(df) <- c("chr", "start", "strand", "type", "pvalue1", "M_s1", "U_s1", "M_s2", "U_s2")
  df$end <- df$start
  df$pvalue2 <- 1 - df$pvalue1
  #filter low coverage sites
  df <- df[which((df$M_s1 + df$U_s1) >= mincov & (df$M_s2 + df$U_s2) >= mincov),]
  df$beta_s1 <- round(df$M_s1 / (df$M_s1 + df$U_s1), 6)
  df$beta_s2 <- round(df$M_s2 / (df$M_s2 + df$U_s2), 6)
  df$FDR1 <- p.adjust(p = df$pvalue1, method = "BH")
  df$FDR2 <- p.adjust(p = df$pvalue2, method = "BH")
  df <- df[,c("chr", "start", "end", "type", "pvalue1", "FDR1", "pvalue2", "FDR2", "M_s1", "U_s1", "M_s2", "U_s2")]
  #get significant subsets
  df1 <- df[which((df$beta_s2 - df$beta_s1 > perc_cut) & (df$FDR1 < FDR)),]
  df2 <- df[which((df$beta_s1 - df$beta_s2 > perc_cut) & (df$FDR2 < FDR)),]
  lst[[samples[1]]] <- df1
  lst[[samples[2]]] <- df2
  return(lst)
}

args = commandArgs(trailingOnly=TRUE)
sample = args[1]
file = args[2]
perc_cut = args[3]
FDR = args[4]
mincov = args[5]

res_lst <- dmcpgs(file, sample, perc_cut, FDR, mincov)
write.table(res_lst[[1]], file = paste0(names(res_lst)[1], ".dmcpgs"), row.names = F, quote = F)
write.table(res_lst[[2]], file = paste0(names(res_lst)[2], ".dmcpgs"), row.names = F, quote = F)

