#' seurat_to_loom
#'
#' Creates a loom function from a Seurat object
#'
#' @param so Seurat object
#' @param filename path to file
#'
#' @return connection to the newly created loom file
#' @export

seurat_to_loom <- function(so, filename) {
  lfile <- SeuratDisk::as.loom(so, filename = filename)
  lfile$close_all()
  lfile <- loomR::connect(filename, skip.validate = TRUE, mode = "r+")
  return(lfile)
}

#' Export metadata from loom for the Single Cell Portal
#'
#' @param so Seurat object
#' @param filename file path for export
#' @param clustering_assay string specifying clustering assay
#' @param species__ontology_label string specifying species__ontology_label required by the Single Cell Portal (default: "Mus musculus")
#' @param species string specifying species required by the Single Cell Portal (default: NA)
#' @param disease__ontology_label etc. (default: "diabetes mellitus")
#' @param disease etc. (default: "MONDO_0005015")
#' @param organ__ontology_label etc. (default: "skin epidermis")
#' @param organ etc. (default: "UBERON_0001003")
#' @param library_preparation_protocol__ontology_label etc. (default: "10x 5' v2")
#' @param library_preparation_protocol etc. (default: "EFO_0009900")
#'
#' @return connection to the newly created loom file
#' @export


seurat_add_to_loom_init <- function(so, filename,
                                    clustering_assay = "RNA",
                                    species__ontology_label = "Mus musculus",
                                    species = NA,
                                    disease__ontology_label = "diabetes mellitus", #skin cancer
                                    disease = "MONDO_0005015", # MONDO_0002898
                                    organ__ontology_label = "skin epidermis",
                                    organ = "UBERON_0001003",
                                    library_preparation_protocol__ontology_label = "10x 5' v2", #10x 3' v3
                                    library_preparation_protocol = "EFO_0009900"
) {

  lfile <- SeuratDisk::as.loom(so, filename = filename)
  lfile$close_all()
  lfile <- loomR::connect(filename, skip.validate = TRUE, mode = "r+")

  if (!("scale.data" %in% names(lfile[["layers"]]))) {
    message("scale.data not found. Adding...")
    lfile$add.layer(list(scale.data = so@assays[[clustering_assay]]@scale.data))
  }

  metadata_df <- so@meta.data %>%
    tibble::as_tibble(rownames = "cell") %>%
    dplyr::select(!!rlang::sym("cell"), !!rlang::sym("orig.ident")) %>%
    dplyr::rename(
      NAME = cell,
      biosample_id = orig.ident
    ) %>%
    dplyr::mutate(
      donor_id = stringr::str_extract(biosample_id, "-[:alpha:]*-"),
      species = ifelse(species__ontology_label == "Mus musculus", "NCBITaxon_10090", species),
      species__ontology_label = species__ontology_label,
      disease = disease,
      disease__ontology_label = disease__ontology_label,
      organ = organ,
      organ__ontology_label = organ__ontology_label,
      library_preparation_protocol = library_preparation_protocol, #EFO_0009922
      library_preparation_protocol__ontology_label = library_preparation_protocol__ontology_label,
      sex = "unknown"
    ) %>%
    dplyr::select(stringr::str_subset(colnames(.), "snn_res"))
  metadata <- metadata_df %>% purrr::map(c)
  lfile$add.col.attribute(metadata)

  seurat_add_to_loom(so, lfile, clustering_assay = clustering_assay)

  return(lfile)
}

#' Export metadata from loom for the Single Cell Portal
#'
#' @param so Seurat object
#' @param lfile connection to loom object
#' @param clustering_assay string specifying clustering assay
#' @param overwrite whether or not to overwrite the loom object (default: FALSE)
#' @param npcs How many PCs to include (default: 30)
#' @param str_subset_regex_cols the string for selecting columns by regular expression (default: "_snn_res")
#'
#' @return nothing, the loom file should be updated
#' @export

seurat_add_to_loom <- function(so, lfile,
                               clustering_assay = "RNA",
                               overwrite = FALSE,
                               npcs = 30,
                               str_subset_regex_cols = "_snn_res") {

  if (so@project.name %>% stringr::str_detect("__")) {
    stop("Error: Seurat project name cannot contain '__' as a substring. Please rename. ")
  }

  # TODO add column to indicate cell inclusion

  gene_df <- tibble::tibble(gene = rownames(so)) %>%
    dplyr::mutate(Selected = gene %in% Seurat::VariableFeatures(so))
  pca_gene_loadings <- so@reductions[[paste0(clustering_assay, "_pca")]]@feature.loadings %>%
    tibble::as_tibble(rownames = "gene") %>%
    dplyr::select(1:31)
  gene_df <- gene_df %>%
    dplyr::left_join(pca_gene_loadings, by = "gene") %>%
    tibble::column_to_rownames("gene")
  colnames(gene_df) <- paste0(so@project.name, "__", colnames(gene_df))

  pca_df <- so@reductions[[paste0(clustering_assay, "_pca")]]@cell.embeddings %>%
    tibble::as_tibble(rownames = "cell") %>% .[1:npcs]
  umap_df <- so@reductions$umap@cell.embeddings %>%
    tibble::as_tibble(rownames = "cell")
  cluster_df <- so@meta.data %>%
    dplyr::select(stringr::str_subset(colnames(.), str_subset_regex_cols)) %>%
    tibble::as_tibble(rownames = "cell")

  cell_df <- tibble::tibble(cell = lfile[["col_attrs/CellID"]][]) %>%
    dplyr::left_join(pca_df, by = "cell") %>%
    dplyr::left_join(umap_df, by = "cell") %>%
    dplyr::left_join(cluster_df, by = "cell") %>%
    tibble::column_to_rownames("cell")
  colnames(cell_df) <- paste0(so@project.name, "__",
                              paste0(clustering_assay, "_pca"), "__",
                              colnames(cell_df))
  if (overwrite) {
    lfile$add.col.attribute(cell_df, overwrite = TRUE)
    lfile$add.row.attribute(gene_df, overwrite = TRUE)
  } else {
    lfile$add.col.attribute(cell_df)
    lfile$add.row.attribute(gene_df)
  }

}



