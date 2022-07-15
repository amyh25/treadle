#' Plot gene violin
#'
#' @param plot_df a data frame with a categorical variable and gene expression
#' @param x categorical variable (symbol)
#' @param gene gene expression variable (symbol)
#'
#' @example
#' plot_gene_violin(plot_df, orig.ident, Cd8a)
#'
#' @return A ggplot object showing violin plot
#' @export
#'

plot_gene_violin <- function(plot_df, x, gene) {
  x <- rlang::enquo(x)
  gene <- rlang::enquo(gene)

  plot_df %>%
    ggplot2::ggplot() +
    ggplot2::aes(!!x, !!gene, fill = !!x) +
    ggplot2::geom_violin(scale = "width")
}

subset_columns <- function(plot_df, subset_cols) {

  if (!is.null(subset_cols)) {
    if (length(names(subset_cols)) == 0) {
      stop("Subset_cols must be a vector named by the relevant columns")
    } else {
      for (i in 1:length(subset_cols)) {
        plot_df <- plot_df %>%
          dplyr::filter(!!rlang::sym(names(subset_cols)[i]) == subset_cols[i])
      }
    }
  }

  return(plot_df)

}

#' Plot gene violin directly from loom file
#'
#' @param lfile a connected loom file
#' @param genes a vector of strings specifying genes to plot
#' @param select_layer a string specifying the
#'  layer to extract data from in the loom file (default: "matrix")
#' @param select_row row in the loom file with gene names (default: "Gene")
#' @param select_cols a string specifying cell-level metadata
#'  to use as a categorical variable. Can be a vector of strings
#'  with additional metadata to attach to the ggplot object
#'  for further manipulation (default: "biosample_id", consistent with the
#'  Single Cell Portal metadata requirements)
#' @param subset_cols columns to subset by. NULL by default
#' @param x_index the index of select columns to use as the x-axis in the
#'  violin plot (default: 1)
#' @param ncol number of columns if multiple genes are specified
#'  (default: length(genes))
#' @param combine whether to combine the plots or not.
#'
#' @return A ggplot object showing violin plot
#' @export
#'

plot_gene_violin_from_loom <- function(lfile, genes,
                                       select_layer = "matrix",
                                       select_row = "Gene",
                                       select_cols = "biosample_id",
                                       subset_cols = NULL,
                                       x_index = 1,
                                       ncol = length(genes),
                                       combine = TRUE) {

  genes <- rlang::set_names(genes)
  gene_df <- get_genes_from_loom(lfile, genes, select_layer, select_row)

  cell_selected_metadata_df <- purrr::map_dfc(select_cols, ~lfile[[paste0("col_attrs/", ..1)]][])
  colnames(cell_selected_metadata_df) <- select_cols

  plot_df <- dplyr::bind_cols(gene_df, cell_selected_metadata_df)
  plot_df <- subset_columns(plot_df, subset_cols)

  if (length(genes) > 1) {

    plots <- genes %>%
      purrr::map(~{
        message(paste0("plotting ", ..1))
        plot_gene_violin(plot_df,
                         !!rlang::sym(select_cols[x_index]),
                         !!rlang::sym(..1))
      })

    if (combine) {
      plots <- cowplot::plot_grid(plotlist = plots, ncol = ncol)
    }

  } else {
    plots <- plot_gene_violin(plot_df,
                              !!rlang::sym(select_cols[x_index]),
                              !!rlang::sym(genes))
  }


  return(plots)
}

#' Plot umap from data frame
#'
#' @param plot_df A data frame or tibble for plotting
#' @param color_str String specifying column name to color the UMAP
#' @param umap1_str String specifying column name for x-axis (default: "UMAP_1")
#' @param umap2_str String specifying column name for x-axis (default: "UMAP_2")
#' @param legend_pt_size Point size for dots in legend (default: 3)
#' @param pt_size Point size for dots in plot (default: 0.1)
#' @param pt_stroke Stroke size for dots in plot (default: 0.5)
#' @param label If true, plots label of color_str on plot
#' @param label_text_size Size of label text (default: 6)
#' @param reorder_points Whether or not to order the points by increasing
#' values by the color variable (default: TRUE)
#' @param drop_na Whether or not to drop NAs per the color aesthetic;
#'  also drops rows with "NA" strings
#'
#' @return A ggplot object showing UMAP plot
#' @export
#'

plot_umap <- function(plot_df,
                      color_str = NULL,
                      umap1_str = "UMAP_1",
                      umap2_str = "UMAP_2",
                      legend_pt_size = 3,
                      pt_size = 0.1,
                      pt_stroke = 0.5,
                      label = FALSE, label_text_size = 6,
                      reorder_points = TRUE,
                      drop_na = TRUE) {
  color_sym <- rlang::sym(color_str)

  if (reorder_points & !is.null(color_str)) {
    plot_df <- plot_df %>%
      dplyr::arrange(!!color_sym)
  }
  if (drop_na) {
    plot_df <- plot_df %>%
      dplyr::filter(!!color_sym != "NA") %>%
      tidyr::drop_na(!!color_sym)
  }

  p <- plot_df %>%
    ggplot2::ggplot() +
    ggplot2::aes(!!rlang::sym(umap1_str), !!rlang::sym(umap2_str)) +
    ggplot2::geom_point(size = pt_size, stroke = pt_stroke) +
    ggplot2::coord_fixed() +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = legend_pt_size)))

  if (!is.null(color_str)) {
    p <- p + ggplot2::aes(color = !!color_sym)
  }
  if (label) {
    p <- p +
      ggrepel::geom_text_repel(
        ggplot2::aes(!!rlang::sym(umap1_str),
                     !!rlang::sym(umap2_str),
                     label = !!color_sym),
        color = "black", size = label_text_size,
        data = ~..1 %>%
          dplyr::group_by(!!color_sym) %>%
          dplyr::summarise(UMAP_1 = mean(!!rlang::sym(umap1_str)),
                           UMAP_2 = mean(!!rlang::sym(umap2_str)),
                           .groups = "drop")
      )
  }
  return(p)
}


detect_umap_cols <- function(lfile, umap_str = "UMAP") {

  message("...Automatically detecting UMAP cols...")
  all_umap_names <- names(lfile$col.attrs) %>%
    stringr::str_subset(umap_str)
  umap_cols <- sort(all_umap_names)[1:2]
  message(paste0("UMAP col detected as ", umap_cols, "\n"))

  return(umap_cols)
}


detect_pca_cols <- function(lfile, pc_str = "RNAPC") {

  message("...Automatically detecting PCA cols...")
  all_pca_names <- names(lfile[["col_attrs"]]) %>%
    stringr::str_subset(pc_str)
  pca_cols <- sort(all_pca_names)[1:2]
  message(paste0("PCA col detected as ", pca_cols, "\n"))

  return(pca_cols)
}


#' Plot umap from loom file with gene expression data
#'
#' @param lfile a connected loom file
#' @param genes a vector of genes to plot
#' @param umap_cols a vector of length 2 specifying the names of the umap columns
#'  umap_cols will be automatically detected if not specified
#' @param umap_df a data frame with the umap columns. Will override umap_col specification
#' @param select_layer layer of the loom file to plot expression data from
#'  (default: "matrix")
#' @param select_row row in the loom file with gene names (default: "Gene")
#' @param select_cols a string specifying cell-level metadata
#'  to use as a categorical variable. Can be a vector of strings
#'  with additional metadata to attach to the ggplot object
#'  for further manipulation (default: "biosample_id", consistent with the
#'  Single Cell Portal metadata requirements)
#' @param subset_cols columns to subset by. NULL by default
#' @param ncol number of columns to plot (default: length(genes))
#' @param label whether to label the color variable (default: FALSE;
#' WARNING, currently there's no check on number of categories, so
#' be careful, this may be prone to crashing if you use it on numerical
#' variables)
#' @param combine whether to combine the plots
#' @param drop_na Whether or not to drop NAs per the color aesthetic;
#'  also drops rows with "NA" strings
#' @param pt_size (default: 0.1)
#' @param pt_stroke (default: 0.5)
#'
#' Also will pass any additional named arguments into `plot_umap`
#'
#' @return A ggplot object showing UMAP plot
#' @export
#'

plot_gene_umap_from_loom <- function(lfile, genes,
                                     umap_cols = NULL,
                                     umap_df = NULL,
                                     select_layer = "matrix",
                                     select_row = "Gene",
                                     select_cols = NULL,
                                     subset_cols = NULL,
                                     ncol = length(genes),
                                     label = FALSE,
                                     combine = TRUE,
                                     drop_na = TRUE,
                                     pt_size = 0.1,
                                     pt_stroke = 0.5) {

  genes <- rlang::set_names(genes)

  message("Retreiving genes from loom...")
  gene_df <- get_genes_from_loom(lfile, genes, select_layer, select_row)

  if (!is.null(umap_df)) {

    message("Using umap_df input, ignoring umap_cols. ")
    umap_cols <- colnames(umap_df)

  } else {

    if (is.null(umap_cols)) {
      message("UMAP cols, not specified, detecting...")
      umap_cols <- detect_umap_cols(lfile)
    }
    if (length(umap_cols) != 2) {
      stop("Error: umap_cols must be a character vector of length 2")
    }
    umap_df <- purrr::map_dfc(umap_cols, ~lfile[[paste0("col_attrs/", ..1)]][])
    colnames(umap_df) <- umap_cols

  }

  if (!is.null(select_cols)) {
    cell_selected_metadata_df <- purrr::map_dfc(select_cols, ~lfile[[paste0("col_attrs/", ..1)]][])
    colnames(cell_selected_metadata_df) <- select_cols
    gene_df <- dplyr::bind_cols(gene_df, cell_selected_metadata_df)
  }

  plot_df <- dplyr::bind_cols(gene_df, umap_df)
  plot_df <- subset_columns(plot_df, subset_cols)

  message("plotting genes....")
  if (length(genes) == 1) {
    plots <- plot_umap(plot_df, genes,
                       umap1_str = umap_cols[1], umap2_str = umap_cols[2],
                       drop_na = drop_na,
                       label = label,
                       pt_size = pt_size,
                       pt_stroke = pt_stroke) +
      ggplot2::scale_color_gradient(low = "grey", high = "darkblue")
  } else {

    plot_list <- genes %>%
      purrr::map(~{
        message(paste0("plotting ", ..1))
        plot_umap(plot_df, ..1,
                  umap1_str = umap_cols[1], umap2_str = umap_cols[2],
                  label = label,
                  pt_size = pt_size,
                  pt_stroke = pt_stroke) +
          ggplot2::scale_color_gradient(low = "grey", high = "darkblue")
      })

    if (combine) {
      message("Combining plots...")
      plots <- cowplot::plot_grid(plotlist = plot_list, ncol = ncol)
    } else {
      plots <- plot_list
    }

    message("Done!")

  }

  return(plots)
}

#' Plot umap from loom file with metadata
#'
#' @param lfile a connected loom file
#' @param var_str a cell level variable to plot
#' @param select_cols a string specifying cell-level metadata
#'  to use as a categorical variable. Can be a vector of strings
#'  with additional metadata to attach to the ggplot object
#'  for further manipulation (default: "biosample_id", consistent with the
#'  Single Cell Portal metadata requirements)
#' @param umap_cols a vector of length 2 specifying the names of the umap columns
#'  umap_cols will be automatically detected if not specified
#' @param umap_df a data frame with the umap columns. Will override umap_col specification
#' @param select_layer layer of the loom file to plot expression data from
#'  (default: "matrix")
#' @param select_row row in the loom file with gene names (default: "Gene")
#' @param ncol number of columns to plot (default: length(genes))
#' @param label whether to label the color variable (default: FALSE;
#' WARNING, currently there's no check on number of categories, so
#' be careful, this may be prone to crashing if you use it on numerical
#' variables)
#' @param combine whether to combine the plots
#' @param drop_na Whether or not to drop NAs per the color aesthetic;
#'  also drops rows with "NA" strings
#' @param pt_size (default: 0.1)
#' @param pt_stroke (default: 0.5)
#' @param factor_vars character vector specifying which variables to set as factors (default: NULL)
#'
#' @example
#' plot_var_umap_from_loom(lfile, "RNA_snn_res.0.1", umap_cols = c("UMAP_1", "UMAP_2"))
#' plot_var_umap_from_loom(lfile, "leiden_0.4", umap_df = umap_df, select_row = "var_names")
#'
#' @return A ggplot object showing UMAP plot
#' @export
#'

plot_var_umap_from_loom <- function(lfile, var_str,
                                    select_cols = NULL,
                                    umap_df = NULL,
                                    umap_cols = NULL,
                                    select_layer = "matrix",
                                    select_row = "Gene",
                                    ncol = length(select_cols),
                                    label = FALSE,
                                    combine = TRUE,
                                    drop_na = TRUE,
                                    pt_size = 0.1,
                                    pt_stroke = 0.5,
                                    factor_vars = NULL) {

  var_df <- purrr::map_dfc(var_str, ~lfile[[paste0("col_attrs/", ..1)]][])
  colnames(var_df) <- var_str

  if (!is.null(umap_df)) {

    message("Using umap_df input, ignoring umap_cols. ")
    umap_cols <- colnames(umap_df)

  } else {

    if (is.null(umap_cols)) {
      message("UMAP cols, not specified, detecting...")
      umap_cols <- detect_umap_cols(lfile)
    }
    if (length(umap_cols) != 2) {
      stop("Error: umap_cols must be a character vector of length 2")
    }
    umap_df <- purrr::map_dfc(umap_cols, ~lfile[[paste0("col_attrs/", ..1)]][])
    colnames(umap_df) <- umap_cols

  }

  if (!is.null(select_cols)) {
    cell_selected_metadata_df <- purrr::map_dfc(select_cols, ~lfile[[paste0("col_attrs/", ..1)]][])
    colnames(cell_selected_metadata_df) <- select_cols
    var_df <- dplyr::bind_cols(var_df, cell_selected_metadata_df)
  }

  plot_df <- dplyr::bind_cols(var_df, umap_df)

  if (!is.null(factor_vars)) {
    for (factor_var_str in factor_vars) {
      plot_df[[factor_var_str]] <- factor(plot_df[[factor_var_str]])
    }
  }

  if (length(var_str) == 1) {
    plots <- plot_umap(plot_df, var_str,
                       umap1_str = umap_cols[1], umap2_str = umap_cols[2],
                       drop_na = drop_na,
                       label = label,
                       pt_size = pt_size,
                       pt_stroke = pt_stroke)
  } else {

    plot_list <- var_str %>%
      purrr::map(~plot_umap(plot_df, ..1,
                     umap1_str = umap_cols[1], umap2_str = umap_cols[2],
                     drop_na = drop_na,
                     label = label,
                     pt_size = pt_size,
                     pt_stroke = pt_stroke))

    if (combine) {
      plots <- cowplot::plot_grid(plotlist = plot_list, ncol = ncol)
    } else {
      plots <- plot_list
    }

  }

  return(plots)
}


