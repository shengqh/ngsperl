source('/home/shengq2/program/ngsperl/lib/scRNA/scRNA_func.r')

library(testthat)

test_that("celltype_to_filename", {
  expect_equal(celltype_to_filename("0: Smooth muscle cells"), "0_Smooth_muscle_cells")
  expect_equal(celltype_to_filename("15: Endothelial cells (aorta)"), "15_Endothelial_cells_aorta_")
})