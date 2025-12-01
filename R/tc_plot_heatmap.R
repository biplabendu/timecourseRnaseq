#' Plot a Heatmap with Hierarchical Row Clustering
#'
#' Generates a heatmap using \code{pheatmap::pheatmap} after performing
#' hierarchical clustering on the rows of the provided data. It calculates a
#' zero-centered color palette automatically.
#'
#' @param data A data frame where the first column contains row identifiers
#'   (e.g., gene symbols) and subsequent columns contain numeric data to be
#'   converted to a matrix and plotted.
#' @param cluster_by A character string indicating the dimension to cluster.
#'   Currently, only "row" is supported. Defaults to "row".
#' @param hclust_method A character string specifying the agglomeration method
#'   to be used in \code{\link[stats]{hclust}}. Defaults to "complete".
#' @param n_clusters Integer. The number of clusters to cut the tree into via
#'   \code{\link[stats]{cutree}}. Defaults to 4.
#' @param scale_break_n Integer. The number of color breaks to generate for the
#'   palette. Defaults to 13.
#' @param title A character string for the main title of the plot. Defaults to
#'   "Heatmap".
#' @param treeheight_row Numeric. The height of the row dendrogram in the
#'   heatmap. Defaults to 150.
#' @param col_annot Data frame or \code{NA}. Column annotations passed to the
#'   \code{annotation_col} argument of \code{pheatmap}. Defaults to \code{NA}.
#' @param cluster_cols Logical. Whether to cluster columns. Defaults to
#'   \code{FALSE}.
#' @param ... Additional arguments passed directly to
#'   \code{\link[pheatmap]{pheatmap}}.
#'
#' @return A \code{pheatmap} object representing the generated heatmap.
#'
#' @export
#'
#' @examples
#' # Create example dataset
#' data <- tibble::tribble(
#'   ~plot, ~Y1928, ~Y1929, ~Y1931, ~Y1932, ~Y1933,
#'   "P1",    1.21,   0.68,     0,       0,      0,
#'   "P2",   -2.67,   3.00,  1.07,    1.16,   1.83,
#'   "P3",    1.64,   2.08,  1.82,    1.93,   1.75,
#'   "P4",    2.77,   1.12,  2.14,    0.93,   1.85,
#'   "P5",    3.98,   1.50,  2.73,    2.12,   2.17,
#'   "P6",    0.93,   6.22,  2.28,    0.73,   5.53,
#'   "P7",    1.31,   1.29,  0.990,   0.92,   1.26,
#' )
#'
#' # Plot heatmap
#' tc_plot_heatmap(data, title = "Example Heatmap")
#'
#' # Plot heatmap with row annotations
#' data_col_annot <- data.frame(plot_type = c("high", "high", "high", "low", "low"))
#' rownames(data_col_annot) <- c("Y1928", "Y1929", "Y1931", "Y1932", "Y1933")
#'
#' tc_plot_heatmap(data, n_clusters = 2, col_annot = data_col_annot)
#'
tc_plot_heatmap <- function(data,
                            cluster_by = "row",
                            hclust_method = "complete",
                            n_clusters = 4,
                            scale_break_n = 13,
                            title = "Heatmap",
                            treeheight_row = 150,
                            col_annot = NA,
                            cluster_cols = FALSE,
                            ...) {

  if (cluster_by == "row") {
    # make a list to save the cluster information
    row_clusters <- list()

    # Set genes as rownames and convert it into a matrix
    data <- as.data.frame(data)
    rownames(data) = data[[1]]
    mat_zscore <- as.matrix(data[-1])


    # Hierarchical clustering the geneset
    hclust_rows <- hclust(
      dist(mat_zscore),
      method = hclust_method
    )

    as.dendrogram(hclust_rows) |>
      plot(horiz = TRUE)

    # Make annotations for the heatmaps
    row_clusters <- data.frame(
      cluster = cutree(
        tree = hclust_rows,
        k = n_clusters
      )
    )

    my_pallete <- generate_color_pallete(
      mat = mat_zscore,
      n_colors = scale_break_n
    )

    # Let's plot!
    pheatmap::pheatmap(
      mat_zscore,
      show_rownames = TRUE,
      show_colnames = FALSE,
      ### ROW and COLUMN CLUSTERS ###
      annotation_col = col_annot,
      cluster_cols = cluster_cols,
      annotation_row = row_clusters,
      cutree_rows = n_clusters,
      ### CELL STYLE ###
      gaps_row = c(10),
      border_color = FALSE,
      ### COLOR SCALE ###
      color = my_pallete$p_colors,
      breaks = my_pallete$p_breaks,
      treeheight_row = treeheight_row,
      main = paste0(
        title,
        "\n",
        " n(rows) = ",
        dim(mat_zscore)[1],
        ", n(cols) = ",
        dim(mat_zscore)[2]
      ),
      ## annotation legend
      annotation_legend = TRUE,
      ## Color scale
      legend = TRUE,
      ...
    )
  } else {
    cli::cli_abort("Param `cluster_by` should be set to 'rows'.")
  }
}
