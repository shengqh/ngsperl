# inputFile<-"./BGN.depth.txt.gz"
# name<-"BGN"
# locus<-"chrX:153494939-153509554"
# width<-3000
# height<-3300
# gtf<-"/data/cqs/references/gencode/GRCh38.p13/gencode.v43.annotation.gtf"

library(ggplot2)
library(data.table)
library(dplyr)
library(ggtranscript)
library(patchwork)

if(! exists("inputFile")){
  args <- commandArgs(TRUE)
  if(length(args) == 0){
    setwd("/scratch/cqs/shengq2/ravi_shah_projects/20220107_rnaseq_discovery_hg38/figures")
    inputFile<-"CD7.depth.txt"
    name<-'CD7'
    locus<-"chr17:82314868-82317602"
    width<-1000
    height<-800
  }else{
    inputFile<-args[1]
    name<-args[2]
    locus<-args[3]
    width<-as.numeric(args[4])
    height<-as.numeric(args[5])
  }
}

cat(inputFile, "\n")
dataForPlot<-fread(file=inputFile, sep="\t", header=T)
dataForPlot$FileName<-gsub("_", " ", dataForPlot$File)

outputFile = paste0(name, ".coverage.png")

g1<-ggplot(dataForPlot, aes(Position, File)) +
  geom_tile(aes(fill = PositionCount))+
  scale_fill_gradientn(colours=rev(heat.colors(100)))+
  ylab("") + xlab(paste0(name, " [", locus, "]"))+
  theme_bw() +
  theme(strip.background=element_blank())

if(gtf != ""){
  cat("GTF file:", gtf, "\n")
  all_annotation=fread(gtf, header=F)
  all_exons <- all_annotation |>
    dplyr::filter(V3 == "exon") |>
    dplyr::mutate(
      gene_name = gsub(".*gene_name \"(.*?)\".*", "\\1", V9),
      transcript_name = gsub(".*transcript_name \"(.*?)\".*", "\\1", V9)) |>
    dplyr::rename(
      seqnames=V1,
      type=V3,
      start=V4,
      end=V5,
      strand=V7) |>
    dplyr::select(seqnames,start,end,strand,type,gene_name,transcript_name)

  gene_exons <- all_exons |>
    dplyr::filter(gene_name == name) 

  if (nrow(gene_exons) > 0){
    g2 = ggplot(gene_exons, aes(
            xstart = start,
            xend = end,
            y = transcript_name
        )) + geom_range() +
        geom_intron(
            data = to_intron(gene_exons, "transcript_name"),
            aes(strand = strand)
        ) +
        theme_bw() + 
        theme( axis.title.y = element_blank() )
    
    eachLineHeight = height / (length(unique(dataForPlot$File)) + 3)
    exonHeight = eachLineHeight * (length(unique(gene_exons$transcript_name)) + 1)

    # Get the x-axis range from g1
    x_range <- range(dataForPlot$Position)

    # Apply the same x-axis range to g2
    g2 <- g2 + xlim(x_range)
    g <- g1 / g2 + plot_layout(heights = c(height, exonHeight))      

    height = height + exonHeight
  }
}else{
  g <- g1
}

ggsave(outputFile, g, height=height, width=width, units="px", dpi=300, bg="white")
