setwd('/nobackup/h_cqs/ramirema/vickers/20251205_methylseq_geneplot')

library(data.table)
#BiocManager::install("ggalt")
library(ggalt)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(patchwork)

plg_plg_vs_unt_plg <- fread("/nobackup/vickers_lab/projects/20251117_14130_DNAMethyl_mm10_dragen/MethylKitDiffAnnovar/result/bPlg_bPlg_vs_unt_bPlg.dmcpgs.annovar.final.tsv", data.table = F)

gtf_file <- "/data/cqs/references/gencode/GRCm38.p6/gencode.vM25.annotation.gtf"
gtf_gr <- rtracklayer::import.gff(con = gtf_file, format = "gtf")

gene="Pik3r2"

gene_methyl <- plg_plg_vs_unt_plg[plg_plg_vs_unt_plg$Gene.refGene == gene, ]

gene_gr <- gtf_gr[gtf_gr$gene_name == gene & gtf_gr$type == "gene"]
xlim=c(start(gene_gr) - 500, end(gene_gr) + 500)

track_df = gene_methyl |>
  dplyr::rename(seqnames=chr, score=meth.diff) |>
  dplyr::mutate(HyperIn=gsub("hyper_in_", "", direction)) |>
  dplyr::select(seqnames, start, end, score, HyperIn)

# Create main methylation plot
g_main = ggplot(data = track_df, aes(x = start, y = score)) +
  geom_lollipop(color="gray", 
                point.size=0) +
  geom_point(aes(color=HyperIn), size=2) +
  theme_classic() +
  xlim(xlim) +
  labs(x = paste0("Chromosome ", unique(track_df$seqnames), " Position (bp)"), 
       y="Methylation Difference %", 
       title=gene) + 
  theme(plot.title = element_text(hjust = 0.5))

# Create gene annotation track
gene_transcripts <- gtf_gr[gtf_gr$gene_name == gene & 
                gtf_gr$type %in% c("exon", "UTR")]

# Keep the transcript with highest exon number
exon_counts <- gene_transcripts %>%
  as.data.frame() %>%
  dplyr::filter(type == "exon") %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarise(n_exons = n()) %>%
  dplyr::arrange(desc(n_exons))

top_transcript_id <- exon_counts$transcript_id[1]

gene_transcripts <- gene_transcripts[gene_transcripts$transcript_id == top_transcript_id]

# Prepare annotation data
anno_df <- as.data.frame(gene_transcripts)

g_gene = ggplot() +
  geom_rect(data = anno_df,
            aes(xmin = start, xmax = end, ymin = 0.3, ymax = 0.7, fill = type)) +
  geom_segment(data = as.data.frame(gene_gr), 
               aes(x = start, xend = end, y = 0.5, yend = 0.5), 
               linewidth = 1, color = "black",
               arrow = arrow(length = unit(0.3, "cm"), 
                           ends = ifelse(as.character(strand(gene_gr)) == "+", "last", "first"))) +
  scale_fill_manual(values = c("exon" = "#4CAF50", "UTR" = "lightblue")) +
  theme_void() +
  xlim(xlim) +
  ylim(0, 1) +
  labs(fill = "Feature")


# Combine plots
g = g_main / g_gene + plot_layout(heights = c(3, 1), guides = "collect") & theme(legend.justification = "left")

ggsave(paste0(gene, "_methylation_plot.png"), plot = g, width = 8, height = 4, dpi = 300, units = "in", bg = "white")

