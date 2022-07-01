#' Export metadata from loom for the Single Cell Portal
#'
#' @param lfile a connected loom file
#' @param output_dir output directory
#' @param filename filename (default: "scp_metadata")
#' @param str_subset_regex_group string denoting regular expression to select
#'  group variables from cell-level metadata
#' @param str_subset_regex_numeric string denoting regular expression to select
#'  numeric variables from cell-level metadata
#' @param fill_metadata boolean specifying whether or not to fill in
#'  arbitrary metadata. USE WITH CAUTION and don't forget about it.
#' @param fill_metadata_cellname Cellname to specify when using fill_metadata option
#'
#' This function will create a 'tmp.txt' file. If you routinely
#' name files that... maybe watch out? Will be fixed eventually.
#'
#' @return Nothing. Check the output directory for a new file.
#' @export
#'

scp_export_metadata_from_loom <- function(lfile, output_dir,
                                          filename = "scp_metadata",
                                          str_subset_regex_group = NA,
                                          str_subset_regex_numeric = NA,
                                          fill_metadata = FALSE,
                                          fill_metadata_cellname = NA) {
  metadata_cols <- c(
    "NAME",
    "biosample_id", "donor_id", "sex",
    "disease", "disease__ontology_label",
    "library_preparation_protocol", "library_preparation_protocol__ontology_label",
    "organ", "organ__ontology_label",
    "species", "species__ontology_label"
  )
  metadata_col_types <- c(
    "TYPE",
    "group", "group", "group",
    "group", "group",
    "group", "group",
    "group", "group",
    "group", "group"
  )
  metadata_col_types_line <- paste0(metadata_col_types, collapse = "\t")

  if (fill_metadata) {
    warning("fill_metadata set to TRUE. Filling in arbitrary metadata. You WILL need to change this later.")

    if (is.na(fill_metadata_cellname)) {
      stop("Must specify fill_metadata_cellname to autofill metadata")
    } else {
      message("Checking NAME IDs are unique...")
      name_ids <- lfile$col.attrs[[fill_metadata_cellname]][]
      is_unique <- length(name_ids) == length(unique(name_ids))
      if (is_unique)
        message("NAME IDs are unique!")
      else
        stop("NAME IDs are not unique, select other column. ")
    }

    metadata_df <- tibble::tibble(NAME = name_ids) %>%
      dplyr::mutate(
        biosample_id = "sample",
        donor_id = "donor",
        sex = "unknown",
        disease__ontology_label = "diabetes mellitus",
        disease = "MONDO_0005015",
        library_preparation_protocol__ontology_label = "10x 5' v2",
        library_preparation_protocol = "EFO_0009900",
        organ__ontology_label = "skin epidermis",
        organ = "UBERON_0001003",
        species__ontology_label = "Mus musculus",
        species = "NCBITaxon_10090"
      )

  } else {

    if (!is.na(str_subset_regex_group)) {
      other_cols_to_select_group <- str_subset(names(lfile$col.attrs), str_subset_regex_group)
      metadata_cols <- c(metadata_cols, other_cols_to_select_group)
      metadata_col_types <- c(metadata_col_types,
                              rep("group", (length(other_cols_to_select_group))))
    }
    if (!is.na(str_subset_regex_numeric)) {
      other_cols_to_select_num <- str_subset(names(lfile$col.attrs), str_subset_regex_numeric)
      metadata_cols <- c(metadata_cols, other_cols_to_select_num)
      metadata_col_types <- c(metadata_col_types,
                              rep("numeric", (length(other_cols_to_select_num))))
    }

    metadata_df <- metadata_cols %>%
      set_names() %>%
      purrr::map_dfc(~lfile$col.attrs[[..1]][])

  }

  colnames(metadata_df) <- metadata_df %>%
    colnames() %>%
    stringr::str_replace_all("\\.", "_")

  write_tsv(metadata_df, file.path(output_dir, "tmp.txt"))
  f_clustering_raw <- file(file.path(output_dir, "tmp.txt"),
                           open = "r")
  lines <- readLines(f_clustering_raw)
  close(f_clustering_raw)
  new_lines <- c(lines[1], metadata_col_types_line, lines[2:length(lines)])
  f_clustering <- file(file.path(output_dir, paste0(filename, ".txt")), open = "w")
  writeLines(new_lines, f_clustering)
  close(f_clustering)
  file.remove(file.path(output_dir, "tmp.txt"))
  message("done!")
}


#' Export clustering from loom for the Single Cell Portal
#'
#' @param lfile a connected loom file
#' @param output_dir output directory
#' @param filename filename (default: "scp_clustering")
#' @param str_subset_regex_group string denoting regular expression to select
#'  group variables from cell-level metadata
#' @param str_subset_regex_numeric string denoting regular expression to select
#'  numeric variables from cell-level metadata
#' @param include_cells vector of cells to export. By default, includes
#'  all cells (default: NULL)
#' @param NAME string specifying the cell IDs to go into the NAME column
#'  (default: NA)
#'
#' @return Nothing. Check the output directory for a new file.
#' @export
#'


scp_export_clustering_from_loom <- function(lfile,
                                            output_dir,
                                            filename = "scp_clustering",
                                            umap_cols = NULL,
                                            str_subset_regex_group = NA,
                                            str_subset_regex_numeric = NA,
                                            include_cells = NULL,
                                            NAME = NA) {

  ### get cell cols
  if (is.na(NAME))
    metadata_cols <- c("NAME")
  else
    metadata_cols <- NAME
  metadata_col_types <- c("TYPE")

  ### add UMAP cols
  if (is.null(umap_cols)) {
    message("UMAP cols, not specified, detecting...")
    umap_cols <- detect_umap_cols(lfile)
  }
  if (length(umap_cols) != 2) {
    stop("Error: umap_cols must be a character vector of length 2")
    print(umap_cols)
  }
  metadata_cols <- c(metadata_cols, umap_cols)
  metadata_col_types <- c(metadata_col_types, "numeric", "numeric")

  ### add other cols
  if (!is.na(str_subset_regex_group)) {
    other_cols_to_select_group <- str_subset(names(lfile$col.attrs), str_subset_regex_group)
    metadata_cols <- c(metadata_cols, other_cols_to_select_group)
    metadata_col_types <- c(metadata_col_types,
                            rep("group", (length(other_cols_to_select_group))))
  }
  if (!is.na(str_subset_regex_numeric)) {
    other_cols_to_select_num <- str_subset(names(lfile$col.attrs), str_subset_regex_numeric)
    metadata_cols <- c(metadata_cols, other_cols_to_select_num)
    metadata_col_types <- c(metadata_col_types,
                            rep("numeric", (length(other_cols_to_select_num))))
  }

  ### get data
  metadata_col_types_line <- paste0(metadata_col_types, collapse = "\t")
  metadata_df <- metadata_cols %>% set_names() %>%
    map_dfc(~lfile$col.attrs[[..1]][])

  ### change name of NAME and umap cols
  tbl_names <- names(metadata_df)
  tbl_names[1] <- "NAME"; tbl_names[2] <- "X"; tbl_names[3] <- "Y"
  names(metadata_df) <- tbl_names

  ### select cells
  if (!is.null(include_cells)) {
    metadata_df <- metadata_df %>%
      filter(NAME %in% include_cells)
  }

  ### write
  write_tsv(metadata_df, file.path(output_dir, "tmp.txt"))
  f_clustering_raw <- file(file.path(output_dir, "tmp.txt"),
                           open = "r")
  lines <- readLines(f_clustering_raw)
  close(f_clustering_raw)
  new_lines <- c(lines[1], metadata_col_types_line, lines[2:length(lines)])
  f_clustering <- file(file.path(output_dir, paste0(filename, ".txt")), open = "w")
  writeLines(new_lines, f_clustering)
  close(f_clustering)
  file.remove(file.path(output_dir, "tmp.txt"))
  message("done!")
}


#' Export matrix from loom for the Single Cell Portal
#'
#' @param lfile a connected loom file
#' @param output_dir output directory
#' @param layer_name the name of the layer to export. The file will be named
#'  accordingly
#' @param gzip whether to use gzip to compress the object (default: TRUE)
#' @param NAME string specifying the cell IDs to select (default: "NAME")
#' @param Gene string specifying the gene IDs to select (default: "Gene")
#'
#' @return Nothing. Check the output directory for a new file.
#' @export
#'


scp_export_matrix_from_loom <- function(lfile,
                                        output_dir,
                                        layer_name = "matrix",
                                        gzip = TRUE,
                                        NAME = "NAME",
                                        Gene = "Gene") {
  if (layer_name == "matrix") {
    matrix_t <- t(lfile[["matrix/"]][,])
  } else {
    matrix_t <- t(lfile[[paste0("layers/", layer_name)]][,])
  }
  rownames(matrix_t) <- lfile$row.attrs[[Gene]][]
  colnames(matrix_t) <- lfile$col.attrs[[NAME]][]
  tbl_t <- as_tibble(matrix_t, rownames = "GENE")
  if (gzip) {
    write_tsv(tbl_t, file.path(output_dir, paste0("scp_layer_", layer_name, ".txt.gz")))
  } else {
    write_tsv(tbl_t, file.path(output_dir, paste0("scp_layer_", layer_name, ".txt")))
  }
  message("done!")
}
