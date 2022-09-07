#' plot_umap
#'
#' This functions plots a UMAP from a data frame
#'
#' @param umap_df a tibble containing umap data
#' @param var_str (optional) variable encoded in color of plot, by default will grab the third column in \code{umap_df}
#' @param umap_cols (optional) vector of strings specifying umap columns to plot, by default will grab the first two columns in \code{umap_df}
#' @param pt_size point size passed to \code{geom_scattermore} or \code{geom_point}
#' @param rasterize rasterize with \code{geom_scattermore} or not with \code{geom_point}
#' @param pixels number of pixels if rasterize is TRUE
#' @param reorder_points whether to order points based on expression
#' @param remove_axis removes axis
#' @param remove_axis_labels changes x and y axis labels to UMAP 1 and UMAP 2
#'
#' @return A ggplot2 object
#' @export

plot_umap <- function(umap_df,
                      var_str = NA,
                      umap_cols = NA,
                      pt_size = 2,
                      rasterize = TRUE,
                      pixels = c(1024, 1024),
                      reorder_points = TRUE,
                      remove_axis = TRUE,
                      rename_axis_labels = TRUE) {
  # if no variable specified, grab the third column
  if (is.na(var_str))
    var_str <- colnames(umap_df)[3]

  var_sym <- rlang::sym(var_str)
  if (reorder_points)
    umap_df <- umap_df %>%
    dplyr::arrange(!!var_sym)

  # set umap cols if not specified
  if (is.na(umap_cols))
    umap_cols <- colnames(umap_df)[1:2]
  else {
    if (!(length(umap_cols) == 2))
      stop("umap_cols must be a character vector of length 2")
  }

  # actually plot
  p <- umap_df %>%
    ggplot2::ggplot() +
    ggplot2::aes_string(umap_cols[1], umap_cols[2]) +
    ggplot2::labs(color = var_str)

  if (rasterize) {
    p <- p +
      scattermore::geom_scattermore(ggplot2::aes(color = get(var_str)),
                                    pointsize = pt_size,
                                    pixels = pixels)
  } else {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = get(var_str)), size = pt_size)
  }

  if (remove_axis)
    p <- p + ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank()
    )

  if (rename_axis_labels)
    p <- p + ggplot2::labs(x = "UMAP dim. 1", y = "UMAP dim. 2")

  return(p)
}
