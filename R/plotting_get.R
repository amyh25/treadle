#' get_metadata
#'
#' Get metadata from a connceted loom object
#'
#' @param lfile a connected loom object
#' @param obs (optional) character vector with specified observations; by default will take all
#'
#' @return A tibble with metadata
#' @export

get_metadata <- function(lfile, obs = NA) {
  expected_dim <- l_full[["matrix"]]$dims[1]
  if (is.na(obs))
    obs <- names(lfile$col.attrs)
  metadata_full <- obs %>%
    unname() %>%
    set_names() %>%
    map_dfc(~ {
      if (lfile[[paste0("col_attrs/", ..1)]]$dims[1] == expected_dim) {
        lfile[[paste0("col_attrs/", ..1)]][]
      }
    })
}

#' get_genes
#'
#' Get genes from a connected loom object
#'
#' @param lfile a connected loom object
#' @param genes a vector of strings specifying genes
#' @param select_layer a string specifying the
#'  layer to extract data from in the loom file
#' @param select_row row in the loom file with gene names
#' @param select_cols optional; set obs names
#'
#' @return A tibble with gene expression data
#' @export

get_genes <- function(lfile,
                      genes,
                      select_layer = "matrix",
                      select_row = "Gene",
                      select_cols = NA) {
  # check if row attribute is present
  row_attr_strs <- names(lfile$row.attrs)
  if (!(select_row %in% row_attr_strs))
    stop(paste0(select_row, " not found among row attributes"))

  # check if genes are present
  all_genes <- lfile$row.attrs[[select_row]][]
  for (g in genes) {
    if (!(g %in% all_genes))
      stop(paste0(g, " not found in ", select_row, "."))
  }

  # get indices
  genes_i <- match(genes, all_genes)

  # get values
  if (select_layer == "matrix") {
    gene_df <- tibble::as_tibble(lfile[[select_layer]][, genes_i])
  } else {
    gene_df <-
      tibble::as_tibble(lfile[[paste0("layers/", select_layer)]][, genes_i])
  }

  # set names
  colnames(gene_df) <- genes
  if (!is.na(select_cols)) {
    select_cols_df <- select_cols %>%
      rlang::set_names() %>%
      purrr::map_dfc(~lfile$col.attrs[[..1]][])
    gene_df <- gene_df %>%
      dplyr::bind_cols(select_cols_df)
  }

  return(gene_df)
}

#' get_umap_df
#'
#' Get umap coordinates from a connected loom object
#'
#' @param lfile a connected loom file
#' @param umap_cols specify names of UMAP columns to select; otherwise will be
#'  automatically detected by \code{\link[treadle::detect_umap_cols()]{detect_umap_cols}}
#' @param obs a vector strings specifying observation variables
#' @param genes a vector of strings specifying genes
#' @param ... parameters for \code{\link[treadle::get_genes()]{get_genes}}
#'
#' @return A tibble with gene expression data
#' @export

get_umap_df <- function(lfile, umap_cols = NA, obs = NA, genes = NA, ...) {
  # get umap cols
  if (is.na(umap_cols))
    umap_cols <- detect_umap_cols(lfile)
  umap_df <- umap_cols %>%
    rlang::set_names() %>%
    purrr::map_dfc(~lfile$col.attrs[[..1]][])

  # check obs cols exist
  obs_names <- names(lfile$col.attrs)
  for (ob in obs)
    if (!(ob %in% obs_names)) stop(paste0(ob, " not found in col.attrs"))

  # get obs cols
  obs_df <- obs %>%
    rlang::set_names() %>%
    purrr::map_dfc(~lfile$col.attrs[[..1]][])
  umap_df <- dplyr::bind_cols(umap_df, obs_df)

  # get genes
  if (!is.na(genes)) {
    genes_df <- get_genes(lfile, genes = genes, ...)
    umap_df <- dplyr::bind_cols(umap_df, genes_df)
  }

  return(umap_df)
}
