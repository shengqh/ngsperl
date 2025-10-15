side_by_side_plot <- function(spatial_so, features, fno, spot_size=1.6, alpha=0.5, omin=0.2, omax=1, dpi=300, with_count=F){
  # features can only be a single string, and should exist in spatial_so@meta.data columns
  if (length(features) != 1) {
    stop("features should be a single string")
  }
  if (!features %in% colnames(spatial_so@meta.data)){
    stop("feature name not in meta data")
  }
  
  # remove _RCTD suffix from features
  features_title <- gsub("_RCTD$", "", features)
  n_col_figure <- 2
  p1 <- SpatialFeaturePlot(object=spatial_so, features='nCount_Spatial.008um',
                           pt.size.factor = spot_size,
                           alpha=c(0, 0),
                           image.alpha=1
  ) + theme(legend.position = "none") + ggtitle('Tissue Image')
  
  if (with_count){
    p1a <- SpatialFeaturePlot(object=spatial_so, features='nCount_Spatial.008um',
                              pt.size.factor = spot_size,
                              alpha=c(0.8, 1),
                              image.alpha=0.1
    ) + theme(legend.position = "right") + ggtitle('nCount')
    n_col_figure <- 3
    p1 <- p1 + p1a
    fig_w <- 12
    fig_h <- 4
  }else{
    fig_w <- 10
    fig_h <- 5
  }
  
  # make the lengend title of p2 as features_title
  p2 <- SpatialFeaturePlot(object=spatial_so, features=features,
                           pt.size.factor = spot_size,
                           alpha=c(omin, omax),  # rescale the feature value to range of 0 to 1, 
                           # cols = brewer.pal(n = 11, name = "Spectral"),
                           # then linearly map to alpha vector
                           # spots with lowest feature value will have alpha=0.2
                           # spots with highest feature value will have alpha=1
                           image.alpha=0.1
  ) + ggtitle(features_title) + theme(legend.position = "right") + labs(fill='RCTD') 
  
  p <- p1 + p2 + plot_layout(guides='collect', ncol=n_col_figure)
  ggsave(filename = fno,plot = p, width=fig_w, height=fig_h, dpi=dpi)
  # return(p)
}

plotting_feature_only <- function(spatial_so, features, fno, spot_size=2, omin=1, omax=1, dpi=500, with_count=F){
  # features can only be a single string, and should exist in spatial_so@meta.data columns
  if (length(features) != 1) {
    stop("features should be a single string")
  }
  if (!features %in% colnames(spatial_so@meta.data)){
    stop("feature name not in meta data")
  }
  
  # remove _RCTD suffix from features
  features_title <- gsub("_RCTD$", "", features)
  n_col_figure <- 1
  # make the lengend title of p2 as features_title
  p2 <- SpatialFeaturePlot(object=spatial_so, features=features,
                           pt.size.factor = 16,
                           alpha=c(0.7, 1),  # rescale the feature value to range of 0 to 1,
                           # then linearly map to alpha vector
                           # spots with lowest feature value will have alpha=0.2
                           # spots with highest feature value will have alpha=1
                           image.alpha=0.2,
                           image.scale = "lowres") + 
    ggtitle(features_title) + 
    theme(legend.position = "right") + 
    labs(fill='RCTD') + 
    # scale_fill_viridis()
    scale_fill_gradient2(low = "orange",mid = "turquoise",high = "purple",midpoint = 0.5)
  
  p = p2 + plot_layout(guides='collect', ncol=1)
  ggsave(filename = fno,plot = p, width=5, height=5, dpi=dpi)
  # return(p)
}

Plot_image_and_counts = function(spatial_so,  fno, spot_size=1.6, omin=0.2, omax=1, dpi=300, with_count=F, lo = "grey50", mi = "yellow3", hi = "red"){
  
  # remove _RCTD suffix from features
  n_col_figure <- 2
  p1 <- SpatialFeaturePlot(object=spatial_so, features='nCount_Spatial.008um',
                           pt.size.factor = spot_size,
                           alpha=c(0, 0),
                           image.alpha=1
  ) + theme(legend.position = "none") + ggtitle('Tissue Image')
  
  if (with_count){
    p1a <- SpatialFeaturePlot(object=spatial_so, features='nCount_Spatial.008um',
                              pt.size.factor = spot_size,
                              alpha=c(1, 1),
                              image.alpha=0.2
    ) + theme(legend.position = "right") + ggtitle('nCount') +
      scale_fill_gradient2(low = lo, mid = mi, high = hi, midpoint = 200)
    n_col_figure <- 2
    p1 <- p1 + p1a
    fig_w <- 12
    fig_h <- 4
  }else{
    fig_w <- 10
    fig_h <- 5
  }
  
  
  p <- p1 + plot_layout(guides='collect', ncol=n_col_figure)
  ggsave(filename = fno,plot = p, width=fig_w, height=fig_h, dpi=dpi)
  # return(p)
}

Plotting_tissue_and_decon = function(spatial_so, features, fno, spot_size=1.6, omin=0.2, omax=1, dpi=300, with_count=F, lo = "grey50", mi = "yellow3", hi = "red"){
  
  if (length(features) != 1) {
    stop("features should be a single string")
  }
  if (!features %in% colnames(spatial_so@meta.data)){
    stop("feature name not in meta data")
  }
  
  
  p1 <- SpatialFeaturePlot(object=spatial_so, features='nCount_Spatial_008um',
                           pt.size.factor = spot_size,
                           alpha=c(0, 0),
                           image.alpha=1
  ) + theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) + 
    ggtitle('Tissue Image')
  
  # remove _RCTD suffix from features
  features_title <- gsub("_RCTD$", "", features)
  p2 <- SpatialFeaturePlot(object=spatial_so, features=features,
                           pt.size.factor = 16,
                           alpha=c(1, 1),  # rescale the feature value to range of 0 to 1,
                           # then linearly map to alpha vector
                           # spots with lowest feature value will have alpha=0.2
                           # spots with highest feature value will have alpha=1
                           image.alpha=0,
                           image.scale = "lowres") + 
    ggtitle(features_title) + 
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) + 
    labs(fill='RCTD') + 
    # scale_fill_viridis()
    scale_fill_gradient2(low = lo, mid = mi, high = hi, midpoint = 0.5)
  
  # make the lengend title of p2 as features_title
  p3 <- SpatialFeaturePlot(object=spatial_so, features=features,
                           pt.size.factor = 16,
                           alpha=c(1, 1),  # rescale the feature value to range of 0 to 1,
                           # then linearly map to alpha vector
                           # spots with lowest feature value will have alpha=0.2
                           # spots with highest feature value will have alpha=1
                           image.alpha=0.4,
                           image.scale = "lowres") + 
    ggtitle(paste0("Tissue + ",features_title)) + 
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5), 
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) + 
    labs(fill='RCTD') + 
    # scale_fill_viridis()
    scale_fill_gradient2(low = lo, mid = mi, high = hi, midpoint = 0.5)
  
  layout <- "
AACC
BBCC
"
  
  p =p1 + p2 + p3 + plot_layout(guides='collect', ncol=2, nrow = 2,design = layout,widths = c(2.5,2.5,5), heights = c(2.5,2.5,5))
  ggsave(filename = fno,plot = p, width=8, height=6, dpi=dpi)
  # return(p)
}

Plotting_all_features_facetted = function(spatial_so, features, fno, spot_size=1.6, omin=0.2, 
                                          omax=1, dpi=300, with_count=F, lo = "grey50", mi = "yellow3", hi = "red", ncols = 4, nrows = 4){
  
  plotList <- lapply(
    features,
    function(key) {
      # Need to assign the plot to a variable
      x <- SpatialFeaturePlot(spatial_so, features = key, image.alpha = 0.3,pt.size.factor = 16) + 
        theme(legend.position = "right",
              legend.direction = "horizontal",
              legend.title.position = "top") +
        scale_fill_gradient2(low = lo, mid = mi, high = hi) +
        labs(fill = gsub("_RCTD","",key))
      x
    }
  )
  p = ggarrange(plotlist = plotList, ncol = ncols, nrow = nrows,legend = 'top')
  
  
  ggsave(filename = fno,plot = p, width=ncols*4, height=nrows*4, dpi=dpi)
  # return(p)
}


plot_double_feature <- function(object, feature1, feature2, fn_lb, top_left = "blue", bottom_right = "orange", bottom_left = "white", top_right = "#FF0000", fn_fig=NULL, dpi=300){
  pw_figure <- paste0(pw_figure_root, fn_lb, "/")
  feature1_alt <- gsub("[^[:alnum:]]+", "_", feature1)
  feature2_alt <- gsub("[^[:alnum:]]+", "_", feature2)
  if (is.null(fn_fig) || length(fn_fig) == 0){
    fn_fig <- paste0(pw_figure, "SpatialFeaturePlotBlend_", feature1, "_vs_", feature2, ".pdf")
  }
  
  plot <- suppressWarnings(SpatialFeaturePlotBlend(
    object = object,
    features = c(feature1, feature2),
    combine = TRUE,
    feature_1_alt_name = feature1_alt,
    feature_2_alt_name = feature2_alt,
    assay = NULL,
    bottom_left = bottom_left,
    bottom_right = bottom_right,
    top_left = top_left,
    top_right = top_right
  ))
  ggsave(fn_fig, plot, width=16, height=5, dpi=dpi)
  
}

SpatialFeaturePlotBlend <- function(object, features, combine = TRUE,
                                    feature_1_alt_name = NULL,
                                    feature_2_alt_name = NULL, assay = NULL,
                                    bottom_left = "#000000",
                                    bottom_right = "#FF0000",
                                    top_left = "#00FF00",
                                    top_right = "#FFFF00",
                                    use_seurat_backend = FALSE,
                                    fp_extra_arguments = list(),
                                    sfp_extra_arguments = list())  {
  # Generate a grid of RGB color values given the requested corner colours.
  gen_color_grid <- function(side_length, bottom_left, bottom_right,
                             top_left, top_right) {
    grad_gen <- function(start, end, n = side_length) {
      colfunc <- colorRampPalette(c(start, end))
      return(colfunc(n))
    }
    # x_y = "x to y"; "bl" = "bottom left", etc
    bl_tl <- grad_gen(bottom_left, bottom_right)
    br_tr <- grad_gen(top_left, top_right)
    l <- lapply(seq_len(length(bl_tl)),
                function(i) {
                  start <- bl_tl[i]
                  end <- br_tr[i]
                  new_grad <- grad_gen(start, end)
                })
    return(t(matrix(unlist(l), ncol = side_length, nrow = side_length)))
  }
  custom_color_SpatialDimPlot <- function(cells_obj, image_name,
                                          new_md_column_name,
                                          colors_per_spot, ...) {
    cells_obj[[new_md_column_name]] <- colors_per_spot
    names(colors_per_spot) <- as.character(colors_per_spot)
    p <- SpatialDimPlot(cells_obj, new_md_column_name,
                        cols = colors_per_spot, images = image_name, ...) +
      ggtitle(new_md_column_name) +
      blend_plot_theme
    return(p)
  }
  extract_colors_from_ggplot <- function(p) {
    built <- ggplot_build(p)$data[[1]]
    if (!is.na(built[1, "fill"])) {
      col_to_use <- "fill"
    } else {
      col_to_use <- "colour"
    }
    return(built[, col_to_use])
  }
  if (length(features) != 2) {
    stop(paste(c("Incorrect number of features. ",
                 "Requires two features, received ",
                 length(features))))
  }
  if (!is.null(assay)) {
    DefaultAssay(object) <- assay
  }
  if (length(fp_extra_arguments) > 0) {
    use_seurat_backend <- TRUE
  }
  blend_plot_theme <- theme(legend.position = "none",
                            plot.title = element_text(hjust = 0.5))
  plot_list_outer <- list()
  for (i in Images(object)) {
    cell_barcodes <- Seurat:::CellsByImage(object, images = i,
                                           unlist = TRUE)
    cells_obj_sub <- object[, cell_barcodes]
    images_sub_list <- list(object[[i]])
    names(images_sub_list) <- i
    cells_obj_sub@images <- images_sub_list
    if (!use_seurat_backend) {
      plot_list <- lapply(features,
                          function(feature) {
                            max_color <- ifelse(feature == features[1],
                                                bottom_right, top_left)
                            SpatialFeaturePlot(object, feature,
                                               images = i,
                                               sfp_extra_arguments) +
                              scale_fill_gradient(low = bottom_left,
                                                  high = max_color) +
                              ggtitle(feature) +
                              blend_plot_theme
                          })
      colors_list <- lapply(plot_list, extract_colors_from_ggplot)
      # Now construct the blended plot
      dat <- FetchData(cells_obj_sub, features)
      side_length <- 100
      col_grid <- gen_color_grid(side_length, bottom_left, bottom_right,
                                 top_left, top_right)
      dat_norm <- apply(dat, 2,
                        function(x) {
                          round((side_length - 1) * x / max(x)) + 1
                        })
      colors_list[[3]] <- sapply(seq_len(nrow(dat_norm)),
                                 function(x) {
                                   col_grid[dat_norm[x, 1],
                                            dat_norm[x, 2]]
                                 })
      legend_grid <- expand.grid(seq(from = min(dat[, features[1]]),
                                     to = max(dat[, features[1]]),
                                     length.out = side_length),
                                 seq(from = min(dat[, features[2]]),
                                     to = max(dat[, features[2]]),
                                     length.out = side_length))
      colnames(legend_grid) <- features
      legend_colors <- c(col_grid)
      legend_grid$color <- legend_colors
      names(legend_colors) <- legend_colors
      legend <- ggplot(legend_grid,
                       aes(x = .data[[features[1]]],
                           y = .data[[features[2]]],
                           color = color)) +
        geom_point(shape = 15, size = 1.9) +
        scale_color_manual(values = legend_colors) +
        coord_cartesian(expand = FALSE) +
        theme(legend.position = "none", aspect.ratio = 1,
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 45,
                                         hjust = 1)) +
        xlab(ifelse(is.null(feature_1_alt_name),
                    features[1], feature_1_alt_name)) +
        ylab(ifelse(is.null(feature_2_alt_name),
                    features[2], feature_2_alt_name))
    } else {
      if (top_right != "#FFFF00") {
        warning(paste("Cannot alter color in top right corner when",
                      "use_seurat_backend is TRUE"))
      }
      vis_reduc <- cells_obj_sub@images[[i]]@coordinates[, c(3, 2)]
      colnames(vis_reduc) <- c("vis_1", "vis_2")
      vis_reduc$vis_2 <- -1 * vis_reduc$vis_2
      vis_reduc_mat <- as.matrix(vis_reduc)
      vis_reduc_obj <- CreateDimReducObject(embeddings = vis_reduc_mat,
                                            key = "vis_")
      cells_obj_sub@reductions$vis <- vis_reduc_obj
      seurat_fp <- do.call(FeaturePlot, c(list(object = cells_obj_sub,
                                               features = features,
                                               reduction = "vis",
                                               blend = TRUE,
                                               cols = c(bottom_left,
                                                        bottom_right,
                                                        top_left),
                                               combine = FALSE),
                                          fp_extra_arguments))
      colors_list <- lapply(seurat_fp[1:3], extract_colors_from_ggplot)
      legend <- seurat_fp[[4]]
    }
    names(colors_list) <- c(features, paste0(features[1], "_", features[2]))
    plot_list <- lapply(names(colors_list),
                        function(x) {
                          do.call(custom_color_SpatialDimPlot,
                                  c(list(cells_obj = cells_obj_sub, i,
                                         x, colors_list[[x]]),
                                    sfp_extra_arguments))
                        })
    plot_list[[4]] <- wrap_plots(ggplot() + theme_void(), legend,
                                 ggplot() + theme_void(), ncol = 1,
                                 heights = c(0.2, 0.6, 0.2))
    plot_list_outer[[i]] <- plot_list
  }
  if (combine == FALSE) {
    return(plot_list_outer)
  } else {
    plot_list_outer <- lapply(plot_list_outer,
                              function(p) {
                                wrap_plots(p, nrow = 1,
                                           widths = c(0.28, 0.28,
                                                      0.28, 0.16))
                              })
    p <- wrap_plots(plot_list_outer, ncol = 1)
    return(p)
  }
}


plot_all <- function(spatial_so, file_lb, pw_figure, dpi=300, spot_size=1.6, with_count=TRUE,
                     col_lo = col_lo, col_mi = col_mi, col_hi = col_hi, plot_individuals = F){
  # cols <- gsub("[^[:alnum:]]+", "_", colnames(spatial_so@meta.data))
  feature_cols <- grep("_RCTD$", colnames(spatial_so@meta.data), value=TRUE)
  
  
  # QC plot
  ViolinPlot <- VlnPlot(spatial_so, features = "nCount_Spatial_008um", raster = FALSE) + NoLegend() + theme(axis.text = element_text(face = "bold", size = 15))
  fno <- paste0(pw_figure, gsub("-", "_", Sys.Date()),"_Raw_Counts_Violin.png")
  ggsave(fno, ViolinPlot, width = 4, height = 5, dpi = 300)
  
  # plot facet
  fno = paste0(pw_figure, "Facetted_allCT_RCTD_prop.png")
  Plotting_all_features_facetted(spatial_so = spatial_so, fno = fno,features = feature_cols,
                                 spot_size = 16,omin = 1,omax = 1,dpi = 1000,with_count = F, 
                                 lo = "grey", mi = "yellow", hi = "red")
  
  # Plot all predicted cell-types
  dominant <- colnames(norm_weights)[apply(norm_weights,1,which.max)]
  spatial_so$dominant <- dominant
  
  fno = paste0(pw_figure, "predicted_celltypes_per_bin.png")
  plot_all_cell_types(spatial_so, fno)
  
  
  fno = paste0(pw_figure, "Image_and_counts.png")
  Plot_image_and_counts(spatial_so, fno, spot_size=spot_size, omin=0.7, omax=1, dpi=dpi, with_count=with_count,
                        lo = "grey", mi = "yellow", hi = "red")
  
  if(plot_individuals){
    for (ifeat in seq_along(feature_cols)) {
      feature_str <- feature_cols[ifeat]
      fno <- paste0(pw_figure, feature_str,".", col_hi[i], ".png")
      Plotting_tissue_and_decon(spatial_so, feature_str, fno, spot_size=spot_size, omin=0.7, omax=1, dpi=dpi, with_count=with_count,
                                lo = "grey", mi = "yellow", hi = "red")
    }
  }
}



plot_all_cell_types <- function(spatial_so, fno){
  
  # all cell types on one figure
  dominant <- colnames(norm_weights)[apply(norm_weights,1,which.max)]
  spatial_so$dominant <- dominant
  
  p1 <- SpatialDimPlot(
    spatial_so,pt.size.factor = 10,
    group.by    = "dominant",
    label       = F
  ) + theme(legend.position = "none")
  p2 <- get_legend(
    SpatialDimPlot(
      spatial_so,pt.size.factor = 10,
      group.by    = "dominant",
      label       = FALSE,    # no labels, just legend
      repel       = FALSE
    ) + theme(
      legend.position   = "right",
      legend.key.width  = unit(0.2, "cm"),
      legend.key.height = unit(0.2, "cm"),
      legend.text       = element_text(size = 15),
      legend.title      = element_text(size = 7),
      # remove all plot panels so only the legend shows
      plot.background   = element_blank(),
      panel.background  = element_blank(),
      axis.text         = element_blank(),
      axis.ticks        = element_blank(),
      axis.title        = element_blank()
    )
  )
  p2 <- ggdraw(p2) + theme(plot.margin=margin(5,5,5,5, "pt"))
  combined <- plot_grid(
    p1, p2,
    ncol       = 2,
    rel_widths = c(8, 4)
  )
  
  # tis = SpatialFeaturePlot(object=spatial_so, features='nCount_Spatial_008um',
  #                          pt.size.factor = 0,
  #                          alpha=c(0, 0),
  #                          image.alpha=1
  # ) + theme(legend.position = "none",
  #           plot.title = element_text(hjust = 0.5),
  #           panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) + 
  #   ggtitle('Tissue Image')
  
  #plot(combined)
  ggsave(filename = fno,plot = combined, width=8, height=5, dpi=1000)
  #ggsave(plot = combined, filename = fno, width=8, height=5)
  # feature1 <- 'PTC_RCTD'
  # feature2 <- 'myCAF_RCTD'
  # plt <- plot_double_feature(spatial_so, feature1, feature2, file_lb)  # alreay saved
  
  
}


t<-function(i){
  cat(colnames(Spatial_SOs[[i]]@meta.data), "\n")
}



#############################################################################
# MORE CUSTOMIZABLE ONES
