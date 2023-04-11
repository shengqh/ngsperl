dmcpgs <- function(file,
                   perc_cut = 0.25,
                   FDR = 0.05){
  require(dplyr)
  lst <- list()
  df <- read.table(file, header = F, sep = "\t")
  colnames(df) <- c("chr", "start", "strand", "type", "pvalue1", "M_s1", "U_s1", "M_s2", "U_s2")
  df$end <- df$start
  df$pvalue2 <- 1 - df$pvalue1
  df$beta_s1 <- round(df$M_s1 / (df$M_s1 + df$U_s1), 6)
  df$beta_s2 <- round(df$M_s2 / (df$M_s2 + df$U_s2), 6)
  df$FDR1 <- p.adjust(p = df$pvalue1, method = "BH")
  df$FDR2 <- p.adjust(p = df$pvalue2, method = "BH")
  df <- df[,c("chr", "start", "end", "type", "pvalue1", "FDR1", "pvalue2", "FDR2", "M_s1", "U_s1", "M_s2", "U_s2")]
  #get significant subsets
  df1 <- df %>%
    dplyr::filter((beta_s2 - beta_s1 > perc_cut) & FDR1 < FDR)
  df2 <- df %>%
    dplyr::filter((beta_s1 - beta_s2 > perc_cut) & FDR2 < FDR)
  lst[[1]] <- df1
  lst[[2]] <- df2
  return(lst)
}

args = commandArgs(trailingOnly=TRUE)
file = args[1]
perc_cut = args[2]
FDR = args[3]

basename <- gsub("^(.+)\\.methdiff", "\\1", file)
res_lst <- dmcpgs(file, perc_cut, FDR)
write.table(res_lst[[1]], file = paste0(basename, "_hypoinSample1.dmcpgs"), row.names = F, quote = F)
write.table(res_lst[[2]], file = paste0(basename, "_hypoinSample2.dmcpgs"), row.names = F, quote = F)

