#' Detect UMAP columns
#'
#' @param lfile link to a loom file
#' @param umap_regex_str string regex to detect UMAP cols using `str_subset`
#'
#' @return a vector of strings of length 2 with names of columns
#' @export

detect_umap_cols <- function(lfile, umap_regex_str = "UMAP|umap") {
  all_umap_names <- names(lfile$col.attrs) %>%
    stringr::str_subset(umap_regex_str)
  message("...Automatically detecting UMAP cols...")
  if (length(all_umap_names) == 0) {
    stop(paste0("No cols detected with regex '", umap_regex_str, "'"))
  }
  umap_cols <- sort(all_umap_names)[1:2]
  message(paste0("UMAP col detected as ", umap_cols, "\n"))

  return(umap_cols)
}
