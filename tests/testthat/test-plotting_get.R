test_that("get_genes returns normalized CCL4 expr exactly like Seurat", {
  lfile <- new_pbmc_loom_connection()
  loom_gene_df <- get_genes(lfile, genes = "CCL4")
  lfile$close_all()

  so <- new_pbmc_so()
  so_gene_df <- tibble::tibble(CCL4 = unname(so@assays$RNA@data["CCL4",]))

  expect_true(identical(loom_gene_df, so_gene_df))
})

test_that("get_genes returns CCL4 counts exactly like Seurat", {
  lfile <- new_pbmc_loom_connection()
  loom_gene_df <- get_genes(lfile, genes = "CCL4", select_layer = "counts")
  lfile$close_all()

  so <- new_pbmc_so()
  so_gene_df <- tibble::tibble(CCL4 = unname(so@assays$RNA@counts["CCL4",]))

  expect_true(identical(loom_gene_df, so_gene_df))
})

test_that("get_genes rejects rows that don't exist", {
  lfile <- new_pbmc_loom_connection()
  expect_error(get_genes(lfile, genes = "CCL4", select_row = "QWERTY"),
               regexp = "QWERTY not found")
  lfile$close_all()
})
