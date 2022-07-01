#' Get genes from loom file
#'
#' @param lfile a connected loom file
#' @param genes a vector of strings specifying genes
#' @param select_layer a string specifying the
#'  layer to extract data from in the loom file (default: "matrix")
#' @param select_row row in the loom file with gene names (default: "Gene")
#'
#' @return A tibble with gene expression data
#' @export
#'

get_genes_from_loom <- function(lfile, genes,
                                select_layer = "matrix",
                                select_row = "Gene") {

  if (select_layer == "matrix") {
    gene_df <- map_dfc(
      genes,
      ~{lfile[[select_layer]][, which(lfile[[paste0("row_attrs/", select_row)]][] == ..1)]}
    )
  } else {
    gene_df <- map_dfc(
      genes,
      ~lfile[[paste0("layers/", select_layer)]][, which(lfile[[paste0("row_attrs/", select_row)]][] == ..1)]
    )
  }

  colnames(gene_df) <- genes
  return(gene_df)
}
