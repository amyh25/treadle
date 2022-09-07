new_pbmc_loom_connection <- function() {
  suppressWarnings({
    lfile <- loomR::connect(test_path("testdata", "processed_mini_10x_pbmc.loom"),
                            skip.validate = TRUE)
  })
  return(lfile)
}

new_pbmc_so <- function() {
  so <- readRDS(test_path("testdata", "processed_mini_10x_pbmc.rds"))
  return(so)
}
