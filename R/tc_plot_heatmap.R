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
#' @param cluster_rows boolean values determining if rows should be clustered
#'   (Default: TRUE)
#' @param gaps_row vector of row indices that show where to put gaps into
#'   heatmap. Set to NA if the rows are not clustered. (Default: 10)
#' @param hclust_method A character string specifying the agglomeration method
#'   to be used in \code{\link[stats]{hclust}}. Defaults to "complete".
#' @param pheatmap_method A character string specifying the agglomeration method
#'   to be used in \code{\link[pheatmap]{pheatmap}}. If \code{NULL},
#'   defaults to \code{hclust_method}.
#' @param dist_method A character string specifying the distance measure to use
#'   clustering rows. Possible values are all distances supported by
#'   \code{\link[stats]{dist}}: "euclidean", "maximum", "manhattan", "canberra",
#'   "binary" or "minkowski". Defaults to "euclidean".
#' @param n_clusters Integer. The number of clusters to cut the tree into via
#'   \code{\link[stats]{cutree}}. Defaults to 4.
#' @param scale_break_n Integer. The number of color breaks to generate for the
#'   palette. Defaults to 13.
#' @param title A character string for the main title of the plot. Defaults to
#'   "Heatmap".
#' @param show_dim_in_title Logical. If true, appends the dimension of dataset
#'    in the title of the plot. Defaults to \code{TRUE}
#' @param treeheight_row Numeric. The height of the row dendrogram in the
#'   heatmap. Defaults to 150.
#' @param col_annot Data frame or \code{NA}. Column annotations passed to the
#'   \code{annotation_col} argument of \code{pheatmap}. Defaults to \code{NA}.
#' @param cluster_cols Logical. Whether to cluster columns. Defaults to
#'   \code{FALSE}.
#' @param plot_dendogram Logical. Whether to print dendogram individually.
#'   Defaults to \code{FALSE}.
#' @param return_cluster_info Logical. If true, returns cluster information as tbl.
#'   Defaults to \code{FALSE}.
#' @param show_rownames boolean specifying if row names are be shown (Default: TRUE).
#' @param show_colnames boolean specifying if column names are be shown (Default: FALSE).
#' @inheritParams generate_color_pallete
#' @param ... Additional arguments passed directly to
#'   \code{\link[pheatmap]{pheatmap}}.
#'
#' @return A \code{tibble} with two columns, one with the unique row identifier,
#' and another specifying the cluster ID assigned to the row identifier.
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
                            cluster_rows = TRUE,
                            gaps_row = c(10),
                            hclust_method = "complete",
                            pheatmap_method = NULL,
                            dist_method = "euclidean",
                            n_clusters = 4,
                            scale_break_n = 13,
                            title = "Heatmap",
                            show_dim_in_title = TRUE,
                            treeheight_row = 150,
                            col_annot = NA,
                            cluster_cols = FALSE,
                            return_cluster_info = TRUE,
                            show_rownames = TRUE,
                            show_colnames = FALSE,
                            color_set = c("maroon", "white", "orange"),
                            plot_dendogram = FALSE,
                            ...) {


  # Set genes as rownames and convert it into a matrix
  data <- as.data.frame(data)
  rownames(data) = data[[1]]
  mat_zscore <- as.matrix(data[-1])


  if (cluster_by == "row") {
    # make a list to save the cluster information
    row_clusters <- list()

    # Hierarchical clustering the geneset
    hclust_rows <- hclust(
      dist(mat_zscore, method = dist_method),
      method = hclust_method
    )

    if (plot_dendogram) {
      as.dendrogram(hclust_rows) |>
        plot(horiz = TRUE)
    }

    # Make annotations for the heatmaps
    row_clusters <- data.frame(
      cluster = cutree(
        tree = hclust_rows,
        k = n_clusters
      )
    )

    my_pallete <- generate_color_pallete(
      mat = mat_zscore,
      color_set = color_set,
      n_colors = scale_break_n
    )

    if (show_dim_in_title) {
      plot_title <- paste0(
        title,
        "\n",
        " n(rows) = ",
        dim(mat_zscore)[1],
        ", n(cols) = ",
        dim(mat_zscore)[2]
      )
    } else {
      plot_title = title
    }

    # Set params if cluster_rows == FALSE
    if (!cluster_rows) {
      gaps_row = NULL
      row_clusters = NA
    }

    # Let's plot!
    p <- pheatmap::pheatmap(
      mat_zscore,
      show_rownames = show_rownames,
      show_colnames = show_colnames,
      ### ROW and COLUMN CLUSTERS ###
      annotation_col = col_annot,
      cluster_cols = cluster_cols,
      cluster_rows = cluster_rows,
      annotation_row = row_clusters,
      cutree_rows = n_clusters,
      ### CELL STYLE ###
      gaps_row = gaps_row,
      border_color = FALSE,
      ### COLOR SCALE ###
      color = my_pallete$p_colors,
      breaks = my_pallete$p_breaks,
      treeheight_row = treeheight_row,
      main = plot_title,
      ## annotation legend
      annotation_legend = TRUE,
      ## Color scale
      legend = TRUE,
      clustering_method = ifelse(
        !is.null(pheatmap_method),
        pheatmap_method,
        hclust_method
      ),
      clustering_distance_rows = dist_method,
      ...
    )

    if (return_cluster_info & cluster_rows) {
      row_clusters |>
        tibble::rownames_to_column("row") |>
        tibble::as_tibble() |>
        dplyr::arrange(
          cluster
        )
    }

  } else {
    cli::cli_abort("Param `cluster_by` should be set to 'rows'.")
  }
}
