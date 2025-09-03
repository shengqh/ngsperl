#Deconvolution functions
savefig <- function(p, filename, width=5, height=5, dpi=300){
  ggsave(filename, p, width=width, height=height, dpi=dpi)
}

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
  p1 <- SpatialFeaturePlot(object=spatial_so, features='nCount_Spatial',
                           pt.size.factor = spot_size,
                           alpha=c(0, 0),
                           image.alpha=1
  ) + theme(legend.position = "none") + ggtitle('Tissue Image')
  
  if (with_count){
    p1a <- SpatialFeaturePlot(object=spatial_so, features='nCount_Spatial',
                              pt.size.factor = spot_size,
                              alpha=c(0.8, 1),
                              image.alpha=0.5
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
                           # then linearly map to alpha vector
                           # spots with lowest feature value will have alpha=0.2
                           # spots with highest feature value will have alpha=1
                           image.alpha=alpha
  ) + ggtitle(features_title) + theme(legend.position = "right") + labs(fill='RCTD')
  
  p <- p1 + p2 + plot_layout(guides='collect', ncol=n_col_figure)
  savefig(p, fno, width=fig_w, height=fig_h, dpi=dpi)
  # return(p)
}

log_msg <- function(message, log_file) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0("[", timestamp, "] ", message, "\n"), file = log_file, append = TRUE)
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

write_seurat_preprocessing = function(samples, groupings, where_to){
  write(paste0(
    "library(Seurat)
     library(dplyr)
     library(spacexr)
     library(Matrix)
     library(tidyverse)
     library(hdf5r) # required to read in data file
     library(ggplot2)
     library(patchwork)
     library(cowplot)
     library(grid)
    
    data_dir <- c(",
    paste0("'",samples, "'", collapse = ","),
    ")
    
    
    names(data_dir) <- c(", paste0("'",groupings,"'", collapse = ","),
    ")
    data_names <- names(data_dir)
    
    Spatial_SOs <- list()
    for(i in seq_along(data_dir)){
      
      # force load the high-resolution image
      # img <- Read10X_Image(paste(dat.dir,'spatial',sep='/'),image.name = 'tissue_hires_image.png')
      # sotmp <- Load10X_Spatial(dat.dir,image=img, slice='slice1')
      
      Spatial_SOs[[i]] <- Seurat::Load10X_Spatial(bin.size = 8, data.dir = data_dir[i], slice = 'slice1') 
      Spatial_SOs[[i]]$orig.ident <- names(data_dir)[[i]]
    }
    
    save(Spatial_SOs, file = '",where_to,"result/Spatial_SOs_preprocessed.RData')"), file = paste0(where_to, "pbs/seurat_preprocessing.R"))
}


plot_all <- function(spatial_so, file_lb, dpi=300, spot_size=1.6, with_count=TRUE){
  # cols <- gsub("[^[:alnum:]]+", "_", colnames(spatial_so@meta.data))
  feature_cols <- grep("_RCTD$", colnames(spatial_so@meta.data), value=TRUE)
  pw_figure <- paste0(pw_figure_root, file_lb, "/")
  
  
  # QC plot
  ViolinPlot <- VlnPlot(spatial_so, features = "nCount_Spatial", raster = FALSE) + NoLegend() + theme(axis.text = element_text(face = "bold", size = 15))
  fno <- paste0(pw_figure, "24-1104_Raw_Counts_Violin.png")
  ggsave(fno, ViolinPlot, width = 4, height = 5, dpi = 300)
  
  # Spatial raw counts plot
  fno <- paste0(pw_figure, "24-1104_Raw_Counts_Spatial.png")
  SpatialPlot <- SpatialFeaturePlot(spatial_so, features = "nCount_Spatial") + theme(legend.position = "right", legend.text = element_text(face = "bold", size = 15), legend.title = element_text(face = "bold", size = 15))
  ggsave(fno, SpatialPlot, width = 7, height = 5, dpi = 300)
  
  
  # plot for each feature
  
  for (ifeat in seq_along(feature_cols)) {
    feature_str <- feature_cols[ifeat]
    fno <- paste0(pw_figure, feature_str, ".png")
    side_by_side_plot(spatial_so, feature_str, fno, spot_size=spot_size, alpha=0.6, omin=0.7, omax=1, dpi=dpi, with_count=with_count)
  }
  # plot for double feature
}


plot_all_cell_types <- function(spatial_so){
  
  # all cell types on one figure
  dominant <- colnames(norm_weights)[apply(norm_weights,1,which.max)]
  spatial_so$dominant <- dominant
  
  p1 <- SpatialDimPlot(
    spatial_so,
    group.by    = "dominant",
    label       = F
  ) + theme(legend.position = "none")
  p2 <- get_legend(
    SpatialDimPlot(
      spatial_so,
      group.by    = "dominant",
      label       = FALSE,    # no labels, just legend
      repel       = FALSE
    ) + theme(
      legend.position   = "right",
      legend.key.width  = unit(0.2, "cm"),
      legend.key.height = unit(0.2, "cm"),
      legend.text       = element_text(size = 5),
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
  savefig(combined, 'cell_type_dimplot_new.png', width=8, height=5)
  
  feature1 <- 'PTC_RCTD'
  feature2 <- 'myCAF_RCTD'
  plt <- plot_double_feature(spatial_so, feature1, feature2, file_lb)  # alreay saved
  
  
}

t<-function(i){
  cat(colnames(Spatial_SOs[[i]]@meta.data), "\n")
}

# A function modified from ggspavis::plotCoords
my_plotCoords<-function (df, x_coord = "x", y_coord = "y", annotate, 
    point_size=0.2, point_shape = 15, 
    legend_title=annotate, legend_position = "right", legend_point_size = 3,
    numeric_feature_colors=rev(hcl.colors(9, "Rocket"))) {
    stopifnot(legend_position %in% c("left", "right", "top", "bottom", "none"))
    p <- ggplot(df, aes(x = get(x_coord), y = get(y_coord), color = get(annotate))) + 
        geom_point(size = point_size, shape = point_shape) + 
        coord_fixed() + 
        theme_bw() + 
        theme(legend.position = legend_position, 
              panel.grid = element_blank(),
              axis.title = element_blank(), 
              axis.text = element_blank(), 
              axis.ticks = element_blank(),
              aspect.ratio=1)

    if (is.numeric(df[[annotate]])) {
        p <- p + scale_color_gradientn(legend_title, colours=numeric_feature_colors) + 
          ggtitle(legend_title) + 
          labs(color = NULL) + 
          theme(plot.title = element_text(hjust = 0.5))
    } else if (is.factor(df[[annotate]]) || is.character(df[[annotate]])) {
        p <- p + 
          guides(color = guide_legend(title=legend_title, override.aes = list(size = legend_point_size))) +
          theme(legend.key.size=unit(0, "lines"))
    }              
    
    p
}
