library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

args=commandArgs(trailingOnly=TRUE)
gsea_file=args[1]
gene_list_file=args[2]
output_prefix=args[3]

if(is.na(output_prefix)){
  output_prefix=gsea_file
}

gsea_data <- fread(gsea_file, data.table=FALSE) |>
  rename(Rank=`RANK IN GENE LIST`,
         ES=`RUNNING ES`) |>
  select(Rank, ES)

symbols=fread(gene_list_file, data.table=FALSE) |>
  select(NAME, SCORE)
symbols$Rank=1:nrow(symbols)
cross_symbols=max(symbols$SCORE[symbols$SCORE < 0])
index=which(symbols$SCORE==cross_symbols)[1]

g2data=gsea_data

if(gsea_data$Rank[1] != 0) {
  gsea_data = rbind(data.frame(Rank=0, ES=0), gsea_data)
}
if(gsea_data$Rank[nrow(gsea_data)] != nrow(symbols)) {
  gsea_data = rbind(gsea_data, data.frame(Rank=nrow(symbols), ES=0))
}

max_x=ceiling(max(gsea_data$Rank) / 1000) * 1000

g1=ggplot(gsea_data, aes(x = Rank, y = ES)) +
  geom_line(color = "green", linewidth = 1) +
  geom_hline(yintercept = 0, color = "black") +
  theme_bw() +
  xlim(0, max_x) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.text.x=element_blank(),
    axis.title=element_blank(),
    axis.ticks.x=element_blank(),
  ) +
  labs(
    title = "Enrichment Plot:\nHALLMARK_INTERFERON_ALPHA_RESPONSE"
  )

lower_limit <- quantile(symbols$SCORE, 0.05)
upper_limit <- quantile(symbols$SCORE, 0.95)

g2=ggplot() +
  geom_segment(data=symbols, aes(x=Rank, y=0, xend=Rank, yend=0.4, color=SCORE), size=0.2) + 
  geom_segment(data=g2data, aes(x=Rank, y=0, xend=Rank, yend=1), color="black", size=0.2) + 
  xlim(0, max_x) +
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0,
                         limits = c(lower_limit, upper_limit),
                         oob = scales::squish) +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

g3=ggplot(symbols) +
  geom_segment(aes(x=Rank, y=0, xend=Rank, yend=SCORE), color="gray", size=0.2) + 
  geom_vline(xintercept = index, color="gray", size=0.5, linetype = "dotted") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  labs(
    x = "Rank in Ordered Dataset",
  ) +
  xlim(0, max_x) +
  scale_x_continuous(breaks = seq(0, max(symbols$Rank), by = 2000)) +
  theme(axis.title.y = element_blank(),
        plot.title = element_blank()) +
  annotate("text", x = 0, y = max(symbols$SCORE), 
           label = "'pos' (positively correlated)", 
           color = "red", size = 5, hjust = 0) +
  annotate("text", x = index, y = 0,
           label = paste0("Zero cross at ", index - 1), 
           color = "black", size = 5, hjust = 0.5, vjust=-1) +
  annotate("text", x = max(symbols$Rank), y = min(symbols$SCORE),  # Adjusted position for right-aligned text
           label = "'neg' (negatively correlated)", 
           color = "blue", size = 5, hjust = 1, vjust=0)

legend_data <- data.frame(
  x = 1:3,
  y = 1,
  label = c("Enrichment profile", "Hits", "Ranking metric scores")
)

# Create a dummy plot
dummy_plot <- ggplot(legend_data) +
  geom_segment(aes(x = 1, xend = 1.5, y = 1, yend = 1, color = label), linewidth = 1.5, show.legend = TRUE) +
  scale_color_manual(
    values = c("Enrichment profile" = "green", "Hits" = "black", "Ranking metric scores" = "gray"),
    name = NULL
  ) +
  theme(
    legend.text=element_text(size=12),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black", size = 1),
    legend.key = element_blank()
  )

# Extract and display the legend
legend <- cowplot::get_plot_component(dummy_plot, 'guide-box-top', return_all = TRUE)

ylabel <- "Ranked list metric (PreRanked)    Enrichment Score (ES)"

combined_plot=plot_grid(g1, g2, g3, nrow=3, rel_heights=c(4,1,3), align = "v")

combined_plot2 <- plot_grid(
  combined_plot, legend, 
  ncol = 1, 
  rel_heights = c(1, 0.1)
)

output_file=paste0(output_prefix, ".pdf")
#png(paste0(output_prefix, ".png"), width=6, height=6, units="in", res=300)
pdf(output_file, width=6, height=6)
ggdraw() +
  draw_plot(combined_plot2, x = 0.03, y = 0, width = 0.97, height = 1) +
  draw_label(ylabel, x = 0.02, y = 0.49, angle = 90)
dev.off()

cat("Output written to ", output_file, "\n")

