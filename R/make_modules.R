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
    min_expression = 10,
    min_timepoints = 8,
    method = "wgcna",
    sim_method = "kendall",
    soft_power = NULL,
    min_module_size = 50,
    merge_cutoff_similarity = 0.9
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
  tmp_data <- data |>
    log2_transform_data(
      id_column = id_column,
      log2
    ) |>
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

    # 2. Calculate Kendall's tau-b correlation for each gene-gene pair
    cat("---------------------------------------------------\n")
    cat("2. Calculate similarity of expression \n")
    cat("---------------------------------------------------\n")
    sim_matrix <- calculate_gene_gene_sim(
      data = datExpr,
      method = sim_method,
      cache = FALSE
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
      soft_power <- sft$powerEstimate |> as.integer()
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
      plot = T
    )
    cat("Done.")
    cat("\n\n")

    # 5. Identify modules or clusters
    cat("---------------------------------------------------\n")
    cat("5. Identify modules (clusters) \n")
    cat("---------------------------------------------------\n")
    modules <- create_modules(
      tree = geneTree,
      dissTOM = dissTOM,
      data = datExpr,
      merge_cutoff_similarity = merge_cutoff_similarity,
      min_module_size = min_module_size
    )

    # Create output
    out <- list(
      tree = geneTree,
      dissTOM = dissTOM,
      modules = modules
    )

    return(out)
  } else {
    cli::cli_abort(
      c(
        "Unknown method provided.",
        "i" = "Only method `wgcna` is currently available."
      )
    )
  }

}
