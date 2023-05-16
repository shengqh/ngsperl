require(GenomicRanges, lib.loc = "/data/cqs/ywang/soft/R_4.1/Packages")
require(tidyverse, lib.loc = "/data/cqs/ywang/soft/R_4.1/Packages")
require(methylKit, lib.loc = "/data/cqs/ywang/soft/R_4.1/Packages")
require(ggpubr, lib.loc = "/data/cqs/ywang/soft/R_4.1/Packages")
require(Cairo, lib.loc = "/data/cqs/ywang/soft/R_4.1/Packages")

params <- read.table(parSampleFile2, sep = "\t", header = F)
params <- params %>% column_to_rownames("V2") %>% t() %>% data.frame()
project <- params$task_name
assembly = params$org
var = params$var
mincov = as.numeric(params$mincov)

#prepare to find all specific format files
input_path <- paste("..", "..", "methylkitprep", "result", sep="/")
cpg.pat <- ".CpG.txt$"

cpg.all <- list.files(path = input_path, pattern = cpg.pat, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
cpg.infile <- as.list(paste(input_path, cpg.all, sep = "/"))
file.id <- gsub("(.*)\\.CpG\\.txt$", "\\1", cpg.all)

meta <- read.table(paste0("../../", project, "_meta.tsv"), sep = "\t", header = T)
meta <- meta[order(row.names(meta)),]
#treatment <- as.numeric(factor(meta[,var])) - 1
treatment <- rep(0, nrow(meta))

cpg.obj <- methRead(cpg.infile,
                    sample.id = as.list(file.id),
                    assembly = assembly,
                    treatment = treatment,
                    context = "CpG",
                    mincov = mincov)

filtered.obj <- filterByCoverage(cpg.obj,
                                 lo.count=mincov,
                                 lo.perc=NULL,
                                 hi.count=NULL,
                                 hi.perc=99.99)
#Merging samples
filtered.cpg.meth <- methylKit::unite(filtered.obj, destrand=FALSE)
saveRDS(filtered.cpg.meth, paste0(project, ".filtered.cpg.meth.rds"))

cpg_bvalue_df <- data.frame(filtered.cpg.meth) %>%
  mutate(id = paste(chr, start, end, sep = "_"))
for(i in 1:length(file.id)){
  id <- file.id[i]
  ncs <- paste0("numCs", i)
  cov <- paste0("coverage", i)
  cpg_bvalue_df[, id] <- cpg_bvalue_df[, ncs] / cpg_bvalue_df[, cov]
}
cpg_bvalue_df <- cpg_bvalue_df %>%
  dplyr::select(one_of("id", file.id)) %>%
  column_to_rownames("id")

cpg_bvalue_cor <- cor(cpg_bvalue_df, method = "pearson")

cpg_bvalue_mds <- (1 - cpg_bvalue_cor) %>%
  cmdscale() %>%
  as_tibble()

colnames(cpg_bvalue_mds) <- c("Dim.1", "Dim.2")
#rownames(cpg_bvalue_mds) <- colnames(cpg_bvalue_df)
cpg_bvalue_mds[,var] <- meta[,var]

# Plot MDS
dms_plot <- ggscatter(cpg_bvalue_mds, x = "Dim.1", y = "Dim.2",
                      label = colnames(cpg_bvalue_df),
                      xlab = "MDS 1",
                      ylab = "MDS 2",
                      color = var,
                      palette = "jco",
                      font.label = 8,
                      size = 1,
                      ellipse = F,
                      ellipse.type = "norm",
                      repel = TRUE) +
  theme_bw()


CairoPDF(file = paste0(project, "_methyl_CpG_bvalue_corr_MDS_plot.pdf"), width = 4, height = 3, pointsize = 8, onefile = T)
dms_plot
dev.off()

CairoPNG(file=paste0(project, "_methyl_CpG_bvalue_corr_MDS_plot.png"), height=1500, width=1500, res=300)
print(dms_plot)
dev.off()
