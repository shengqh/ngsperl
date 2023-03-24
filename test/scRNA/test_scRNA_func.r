source('/home/shengq2/program/ngsperl/lib/scRNA/scRNA_func.r')

library(testthat)

test_that("celltype_to_filename", {
  expect_equal(celltype_to_filename("0: Smooth muscle cells"), "0_Smooth_muscle_cells")
  expect_equal(celltype_to_filename("15: Endothelial cells (aorta)"), "15_Endothelial_cells_aorta_")
})

test_that("read_bubble_genes_1", {
  bubblemap_file = "/home/shengq2/program/collaborations/alexander_gelbard/Gelbard_BubbleMap_20220804.xlsx"
  genes_df = read_bubble_genes(bubblemap_file)
  expect_equal(nrow(genes_df), 57)
  expect_equal(genes_df$gene[3], "CD19")
  expect_equal(as.character(genes_df$cell_type[3]), "Pan B Cells")
})

test_that("read_bubble_genes_2", {
  bubblemap_file="/home/shengq2/program/projects/justin_turner/20230107_bubble_airway.xlsx"
  genes_df = read_bubble_genes(bubblemap_file)
  expect_equal(nrow(genes_df), 68)
  expect_equal(genes_df$gene[3], "MS4A1")
  expect_equal(as.character(genes_df$cell_type[3]), "Pan B Cells")
})

test_that("read_bubble_genes_3", {
  bubblemap_file="/home/shengq2/program/collaborations/scRNA/20220707_BubbleMap_mouse_biolegend_aorta.xlsx"
  genes_df = read_bubble_genes(bubblemap_file, species="Mm")
  expect_equal(nrow(genes_df), 36)
  expect_equal(genes_df$gene[1], "Cd68")
  expect_equal(as.character(genes_df$cell_type[1]), "Macrophages")
  expect_equal(genes_df$gene[2], "Adgre1")
  expect_equal(as.character(genes_df$cell_type[2]), "Macrophages")
  expect_equal(genes_df$gene[3], "Cd11b")
  expect_equal(as.character(genes_df$cell_type[3]), "Monocytes")
})


