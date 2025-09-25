rm(list=ls()) 
sample_name='A7372_T'
outFile='A7372_T'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/jennifer_pietenpol_projects/20250626_halo/20250918_normal_tumor_target/LOH_3_facets/result/A7372_T')

### Parameter setting end ###

# modified from https://github.com/mskcc/facets-suite/blob/master/run-facets-wrapper.R

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(facets)
    library(facetsSuite)
    library(dplyr)
    library(purrr)
    library(tibble)
    library(ggplot2)
    library(patchwork)
})

myoptions_tbl = read.table(parSampleFile2, header=FALSE, sep="\t", stringsAsFactors=FALSE)
myoptions = split(myoptions_tbl$V1, myoptions_tbl$V2)

args=list()
args$counts_file = read.table(parSampleFile1, header=FALSE)$V1[1]
args$sample_id = sample_name
args$directory = "."
args$everything = TRUE
args$genome = myoptions$genome
args$cval = 500
args$purity_cval = 1000
args$min_nhet = 15
args$purity_min_nhet = 15
args$snp_window_size = 250
args$normal_depth = 35
args$diplogr = NULL
args$seed = 20250923

# Helper functions ------------------------------------------------------------------------------------------------

# Write out
write = function(input, output) {
    write.table(input, file = output, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
}

# Print run details
print_run_details = function(outfile,
                             cval,
                             min_nhet,
                             purity,
                             ploidy,
                             dipLogR,
                             flags = NULL,
                             ...) {
    
    params = c(...)
    
    run_details = data.frame(
        'sample' = sample_id,
        'purity' = signif(purity, 2),
        'ploidy' = signif(ploidy, 2),
        'dipLogR' = signif(dipLogR, 2),
        'facets_version' = as.character(packageVersion('facets')),
        'cval' = cval,
        'snp_nbhd' = args$snp_window_size,
        'min_nhet' = min_nhet,
        'ndepth' = args$normal_depth,
        'genome' = args$genome,
        'seed' = args$seed,
        'flags' = unlist(flags),
        'input_file' = basename(args$counts_file))
    
    if (length(params) > 0) {
        run_details = data.frame(run_details,
                                 'genome_doubled' = params$genome_doubled,
                                 'fraction_cna' = signif(as.numeric(params$fraction_cna), 2),
                                 'hypoploid' = params$hypoploid,
                                 'fraction_loh' = signif(as.numeric(params$fraction_loh), 2),
                                 'lst' = params$lst,
                                 'ntai' = params$ntelomeric_ai,
                                 'hrd_loh' = params$hrd_loh)
    }
    
    write(run_details, outfile)
}

cf_plot=function(facets_output, method = c("em", "cncf"), plotX = FALSE, 
    genome = c("hg19", "hg18", "hg38"), return_object = FALSE) {
    genome = match.arg(genome, c("hg19", "hg18", "hg38"), several.ok = FALSE)
    method = match.arg(method, c("em", "cncf"), several.ok = FALSE)
    snps = facets_output$snps
    segs = facets_output$segs
    if (!plotX) {
        snps = subset(snps, chrom < 23)
        segs = subset(segs, chrom < 23)
    }
    snps = facetsSuite:::get_cum_chr_maploc(snps, genome)
    mid = snps$mid[names(snps$mid) %in% snps$snps$chrom]
    centromeres = snps$centromeres
    snps = snps$snps
    if (method == "em") {
        cols = c((grDevices::colorRampPalette(c("white", "steelblue")))(10), 
            "papayawhip")[round(10 * segs$cf.em + 0.501)]
        my_ylab = "CF (EM)"
    }
    else if (method == "cncf") {
        my_ylab = "CF (CNCF)"
        cols = c((grDevices::colorRampPalette(c("white", "steelblue")))(10), 
            "papayawhip")[round(10 * segs$cf + 0.501)]
    }
    starts = cumsum(c(1, segs$num.mark))[seq_along(segs$num.mark)]
    ends = cumsum(c(segs$num.mark))
    my_starts = snps[starts, "chr_maploc"]
    my_ends = snps[ends, "chr_maploc"]

    cf_df=data.frame(my_starts = unlist(my_starts[,1]), my_ends= unlist(my_ends[,1]))

    cf = ggplot(segs) + 
        geom_rect(data=cf_df, aes(xmin = my_starts, xmax = my_ends, 
                  ymax = 1, ymin = 0), fill = cols, col = "white", size = 0) + 
        scale_x_continuous(breaks = mid, labels = names(mid), 
            expand = c(0.01, 0)) + scale_y_continuous(expand = c(0, 
        0)) + labs(x = NULL, y = my_ylab) + theme_bw() + theme(axis.text.x = element_text(angle = 0, 
        size = 8, color = "black"), axis.text.y = element_text(angle = 0, 
        size = 8, color = "white"), axis.ticks.y = element_line(color = "white"), 
        text = element_text(size = 10), panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), plot.margin = unit(c(0, 
            1, 0.5, 0), "lines"))
    if (return_object == TRUE) {
        cf
    }
    else {
        suppressMessages(print(cf))
    }
}

# Default set of output plots
print_plots = function(outfile,
                       facets_output,
                       cval,
                       genome) {
    
    plot_title = paste0(sample_id,
                        ' | cval=', cval,
                        ' | purity=', round(facets_output$purity, 2),
                        ' | ploidy=', round(facets_output$ploidy, 2),
                        ' | dipLogR=', round(facets_output$dipLogR, 2))
    
    g1=cnlr_plot(facets_output, genome=genome, return_object = TRUE) + ylab("Copy number\nlog ratio")
    #ggsave(paste0(outfile, ".1.png"), g1, width = 850, height = 999, units="px", dpi = 300, bg="white")

    g2=valor_plot(facets_output, genome=genome, return_object = TRUE) + ylab("Variant allele\nlog odds ratio")
    #ggsave(paste0(outfile, ".2.png"), g2, width = 850, height = 999, units="px", dpi = 300, bg="white")

    g3=icn_plot(facets_output, method = 'em', genome=genome, return_object = TRUE) + ylab("Integer copy number\n(EM)")
    #ggsave(paste0(outfile, ".3.png"), g3, width = 850, height = 999, units="px", dpi = 300, bg="white")

    g4=cf_plot(facets_output, method = 'em', genome=genome, return_object = TRUE) + ylab("CF\n(EM)")
    #ggsave(paste0(outfile, ".4.png"), g4, width = 850, height = 999, units="px", dpi = 300, bg="white")
    
    g5=icn_plot(facets_output, method = 'cncf', genome=genome, return_object = TRUE) + ylab("Integer copy number\n(CNCF)")
    #ggsave(paste0(outfile, ".5.png"), g5, width = 850, height = 999, units="px", dpi = 300, bg="white")
    
    g6=cf_plot(facets_output, method = 'cncf', genome=genome, return_object = TRUE) + ylab("CF\n(CNCF)")
    #ggsave(paste0(outfile, ".6.png"), g6, width = 850, height = 999, units="px", dpi = 300, bg="white")

    g=wrap_plots(list(g1, g2, g3, g4, g5, g6), ncol=1, nrow=6, heights = c(1, 1, 1, .15, 1, .15)) + 
        plot_annotation(title = plot_title) & 
        theme(plot.title = element_text(hjust = 0.5, size = 14))
    ggsave(outfile, g, width = 10, height = 8, units="in", dpi = 300, bg="white")
}

# Print IGV-style .seg file
print_igv = function(outfile,
                     facets_output) {
    
    ii = format_igv_seg(facets_output = facets_output,
                        sample_id = sample_id,
                        normalize = T)
    
    write(ii, outfile)
}

# Define facets iteration
# Given a set of parameters, do:
# 1. Run facets
# 2. Generate and save plots
# 3. Print run iformation, IGV-style seg file, segmentation data
facets_iteration = function(name_prefix, ...) {
    params = list(...)
    
    output = run_facets(read_counts = read_counts,
                        cval = params$cval,
                        dipLogR = params$diplogr,
                        ndepth = params$ndepth,
                        snp_nbhd = params$snp_nbhd,
                        min_nhet = params$min_nhet,
                        genome = params$genome,
                        seed = params$seed)
    
    print_igv(outfile = paste0(name_prefix, '.seg'),
              facets_output = output)
    
    print_plots(outfile = paste0(name_prefix, '.png'),
                facets_output = output,
                cval = params$cval,
                genome = params$genome)
    
    output
}

# Run -------------------------------------------------------------------------------------------------------------

# Name files and create output directory
sample_id = args$sample_id
directory = args$directory

# Read SNP counts file
message(paste('Reading', args$counts_file))
read_counts = read_snp_matrix(args$counts_file)
message(paste('Writing to', directory))

name = paste0(directory, '/', sample_id)

purity_output = facets_iteration(name_prefix = paste0(name, '_purity'), 
                                  sample_id = sample_id,
                                  diplogr = args$diplogr,
                                  cval = args$purity_cval,
                                  ndepth = args$normal_depth,
                                  snp_nbhd = args$snp_window_size,
                                  min_nhet = args$purity_min_nhet,
                                  genome = args$genome,
                                  seed = args$seed)

hisens_output = facets_iteration(name_prefix = paste0(name, '_hisens'),
                                  sample_id = sample_id,
                                  diplogr = purity_output$diplogr,
                                  cval = args$cval,
                                  ndepth = args$normal_depth,
                                  snp_nbhd = args$snp_window_size,
                                  min_nhet = args$purity_min_nhet,
                                  genome = args$genome,
                                  seed = args$seed)
                                  
metadata = c(
    map_dfr(list(purity_output, hisens_output), function(x) { arm_level_changes(x$segs, x$ploidy, args$genome)[-5] }),
    map_dfr(list(purity_output, hisens_output), function(x) calculate_lst(x$segs, x$ploidy, args$genome)),
    map_dfr(list(purity_output, hisens_output), function(x) calculate_ntai(x$segs, x$ploidy, args$genome)),
    map_dfr(list(purity_output, hisens_output), function(x) calculate_hrdloh(x$segs, x$ploidy)),
    map_dfr(list(purity_output, hisens_output), function(x) calculate_loh(x$segs, x$snps, args$genome))
)

qc = map_dfr(list(purity_output, hisens_output), function(x) check_fit(x, genome = args$genome)) %>% 
    add_column(sample = sample_id,
                cval = c(args$purity_cval, args$cval), .before = 1)
# Write QC
write(qc, paste0(name, '.qc.txt'))

# Write gene level // use hisensitivity run
gene_level = gene_level_changes(hisens_output, args$genome) %>% 
    add_column(sample = sample_id, .before = 1)
write(gene_level, paste0(name, '.gene_level.txt'))

# Write arm level // use purity run
arm_level = arm_level_changes(purity_output$segs, purity_output$ploidy, args$genome) %>% 
    pluck('full_output') %>% 
    add_column(sample = sample_id, .before = 1)
write(arm_level, paste0(name, '.arm_level.txt'))

print_run_details(outfile = paste0(name, '.txt'),
                  cval = c(args$purity_cval, args$cval),
                  min_nhet = c(args$purity_min_nhet, args$min_nhet),
                  purity = c(purity_output$purity, hisens_output$purity),
                  ploidy = c(purity_output$ploidy, hisens_output$ploidy),
                  dipLogR = c(purity_output$dipLogR, hisens_output$dipLogR),
                  flags = map(list(purity_output$flags, hisens_output$flags), function(x) paste0(x, collapse = '; ')),
                  metadata)

# Write RDS
saveRDS(purity_output, paste0(name, '_purity.rds'))
saveRDS(hisens_output, paste0(name, '_hisens.rds'))
