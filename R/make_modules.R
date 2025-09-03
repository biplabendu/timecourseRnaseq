#' Make modules from timeseries data
#'
#' @description
#' short desc here
#'
#' @param xxx description
#'
#' @export
#'
#'
make_modules <- function(
    data,
    log2 = TRUE,
    id_column = "gene_name",
    min_expression = NULL,
    min_timepoints = NULL,
    method = "wgcna",
    qc = FALSE,
    sim_method = "kendall",
    soft_power = NULL,
    min_module_size = 50,
    max_modules = 15,
    merge_cutoff_similarity = 0.9,
    plot_network = TRUE,
    plot_network_min_edge = 0.65,
    tidy_modules = TRUE,
    cache = FALSE
) {
  chk::chk_data(data)
  # other checks on the structure of the dataset here
  # e.g., `check_data_structure(data)`
  # 1. one row per gene; no duplicates
  # 2. ...

  # A. Data prep -----------------------------
  cat("---------------------------------------------------\n")
  cat("1. Log2-transform and subset \n")
  cat("---------------------------------------------------\n")
  # log2-transform
  if (log2) {
    data <- data |>
      log2_transform_data(
        id_column = id_column,
        log2
      )
  }
  # Estimate defaults
  if (is.null(min_expression)) {
    min_expression <- estimate_min_expression(data, id_column)
    cat("Estimated min_expression =", min_expression, "\n")
  }
  if (is.null(min_timepoints)) {
    min_timepoints <- ceiling( (ncol(data) - 1) * (2/3) )
    cat("Estimated min_timepoints =", min_timepoints, "\n")
  }
  tmp_data <- data |>
    subset_data(
      min_expression,
      min_timepoints,
      id_column = id_column
    )
  cat("\n\n")

  # B. WGCNA -----------------------------------
  if (tolower(method) == "wgcna") {
    # transpose the dataset
    datExpr <- tmp_data |>
      transpose_data()

    # QC: Perform quality check on the dataset
    if (qc == TRUE) {
      datExpr %>%
        check_sample_quality()
      datExpr %>%
        plot_sample_expression()
    }

    # 2. Calculate Kendall's tau-b correlation for each gene-gene pair
    cat("---------------------------------------------------\n")
    cat("2. Calculate similarity of expression \n")
    cat("---------------------------------------------------\n")
    sim_matrix <- calculate_gene_gene_sim(
      data = datExpr,
      method = sim_method,
      cache = cache
    )
    cat("\n\n")

    cat("---------------------------------------------------\n")
    cat("3. Create adjacency matrix \n")
    cat("---------------------------------------------------\n")
    ### Assess and specify the soft-thresholding-power
    sft <- analyze_network_topology(
      data = datExpr,
      plot = TRUE
    )
    if (is.null(soft_power)) {
      soft_power <- estimate_soft_power(sft)
    } else {
      chk::chk_integer(soft_power)
      chk::chk_gt(soft_power, 1)
    }
    cat(
      glue::glue(
        "Setting soft-thresholding power to: { soft_power }."
      )
    )
    cat("\n")

    ### Create the signed adjacency matrix
    cat("Power-transforming the gene-gene similarity matrix...")

    adj_matrix <- WGCNA::adjacency.fromSimilarity(
      sim_matrix,
      power = soft_power,
      type = "signed"
    ) |>
      as.matrix()
    cat("Done.")
    cat("\n\n")

    cat("---------------------------------------------------\n")
    cat("4. Convert into topological overlap matrix (dissTOM) \n")
    cat("---------------------------------------------------\n")
    cat("Creating dissTOM...")
    dissTOM = 1 - WGCNA::TOMsimilarity(
      adj_matrix,
      verbose = 0
    )
    cat("Done.")
    cat("\n")

    ### Call the hierarchical clustering function
    cat("Performing hierarchical clustering on dissTOM...")
    geneTree = perform_hclust(
      data = dissTOM,
      plot = FALSE
    )
    cat("Done.")
    cat("\n\n")

    # 5. Identify modules or clusters
    cat("---------------------------------------------------\n")
    cat("5. Identify modules (clusters) \n")
    cat("---------------------------------------------------\n")
    modules <- create_modules_auto(
      tree = geneTree,
      dissTOM = dissTOM,
      data = datExpr,
      merge_cutoff_similarity = merge_cutoff_similarity,
      min_module_size = min_module_size,
      max_modules = max_modules
    )

    adj_matrix_ME <- calculate_module_module_sim(
      merged_modules = modules[["modules"]]
    )

    # 6. Tidy modules
    cat("---------------------------------------------------\n")
    cat("6. Tidy modules (clusters) \n")
    cat("---------------------------------------------------\n")

    if (plot_network) {
      plot_adj_as_network(
        layout = igraph::layout_as_tree,
        matrix = adj_matrix_ME$ME,
        min_edge = plot_network_min_edge,
        node_label_size = 1.2,
        node_size = 35,
        edge_size = 5,
        node_frame_col = "grey20",
        node_fill_col = "grey80",
        vertex.frame.width = 4
      )
    }

    if (tidy_modules) {
      module_genes <- tidy_modules(
        merged_modules = modules[["colors"]],
        mapping_tbl = adj_matrix_ME$mapping_tbl,
        data = datExpr
      )

      list(
        modules = modules,
        module_genes = module_genes,
        adj_matrix = adj_matrix
      )
    } else {
      list(
        modules = modules,
        adj_matrix = adj_matrix
      )
    }



  } else {
    cli::cli_abort(
      c(
        "Unknown method provided.",
        "i" = "Only method `wgcna` is currently available."
      )
    )
  }

}
